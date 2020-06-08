#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/microfacet.h>
#include <mitsuba/render/fresnel.h>

NAMESPACE_BEGIN(mitsuba)
template <typename Float, typename Spectrum>
class DisneyBRDF final : public BSDF<Float, Spectrum> {
public:
  MTS_IMPORT_BASE(BSDF, m_flags, m_components)
  MTS_IMPORT_TYPES(Texture, MicrofacetDistribution)

  DisneyBRDF(const Properties &props) : Base(props) {
	m_reflectance = props.texture<Texture>("reflectance", .5f);
	m_metallic = props.texture<Texture>("m_metallic", 0.0f);
	m_roughness = props.texture<Texture>("m_roughness", 1.0f);
	m_ior = props.float_("m_ior", 1.0f);
	m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
	m_components.push_back(m_flags);
  }
  Float blender_schlick_fresnel(Float u) const{
	Float m = clamp(1.0f - u, 0.0f, 1.0f);
	Float m2 = m * m;
	return m2 * m2 * m;
  }
  Float blender_calculate_disney_brdf_diffuse_mul(const SurfaceInteraction3f &si, const Float &roughness, BSDFSample3f &bs, Mask &active) const{
	Vector3f world_L = si.sh_frame.to_world(bs.wo);
	Float NgdotL = dot(si.n, world_L);
	active &= NgdotL > 0.f;
  	Float NdotL = Frame3f::cos_theta(bs.wo); //outcome
	active &= NdotL > 0.f;
	Float NdotV = Frame3f::cos_theta(si.wi); //income
	active &= NdotV > 0.f;
	// if (unlikely(none_or<false>(active) ||
	// 			 !ctx.is_enabled(BSDFFlags::DiffuseReflection))){
	//   return 0.f;
	// }
	Vector3f H = normalize(si.wi + bs.wo);
	Float LdotH = dot(bs.wo, H);
	Float FL = blender_schlick_fresnel(NdotL);
	Float FV = blender_schlick_fresnel(NdotV);
	Float Fd90 = 0.5f + 2.0f * LdotH * LdotH * roughness;
	//lerp(a,b,t) b * t + (1-t)*a
	Float Fd = (1.0f * (1.0f - FL) + Fd90 * FL) * (1.0f * (1.0f - FV) + Fd90 * FV);
	Float value = 1.0f / Float(M_PI) * NdotL * Fd;
	return select(active, value, 0.0f);
  }
  Float pbrt_schlickR0FromEta(Float eta) const{
	return sqr(eta - 1) / sqr(eta + 1); //FIXME: check sqr or sqrt
  }
  Spectrum blender_fresnel(const Spectrum & cspec0, const Float &ior, const Vector3f &m, const Vector3f &wo) const{ //check active condition
	Float F0;
	std::tie(F0, std::ignore, std::ignore, std::ignore)= fresnel(Float(1.0f), ior); //ior might be wrong
	Float F0_norm = 1.0f / (1.0f - F0);
	Float FH_f;
	std::tie(FH_f, std::ignore, std::ignore, std::ignore)= fresnel(dot(wo, m), ior); //ior might be wrong
	Float FH = FH_f * F0_norm;
	return lerp(cspec0, Spectrum(1.0f), FH);
  }
  std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
										   const SurfaceInteraction3f &si,
										   Float sample1,
										   const Point2f &sample2,
										   Mask active) const override {
	MTS_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
	if constexpr (is_polarized_v<Spectrum>) {
        Throw("DisneyBRDF: I don't know how to support polarized mode");
    }
	//We might need to use packet!!!
	Spectrum albedo = m_reflectance->eval(si, active);
	Float roughness = m_roughness->eval_1(si, active);
	Float metallic = m_metallic->eval_1(si, active);
	//Float ior = m_ior->eval_1(si, active);
	Float ior = m_ior;
	Float f = max(ior, 1e-5f);
	Float diffuse_weight = (1.0f - clamp(metallic, 0.0f, 1.0f)); //* (1.0 - calmp(transmissionm 0.0, 1.0))
	//Float final_transmission = clamp(transmission, 0.0, 1.0) * (1.0 - clamp(metallic, 0.0, 1.0));
	Float specular_weight = 1.0;// - final_transmission
	Float specular_tint = 0.0f;
	//Float specular = 0.5f; //blender implementation
	
	Float cdlum = luminance(albedo);
	Spectrum ctint = select(cdlum > 0.0f, albedo/cdlum, Spectrum(0.0f));
	//T=Tagent????
	///m_cdlum???
	
	
	Spectrum BSDF(0.0f);
	BSDFSample3f bs = zero<BSDFSample3f>();
	
	bs.wo = warp::square_to_cosine_hemisphere(sample2);

	if (any_or<true>(diffuse_weight > 1e-5f)) {
	  //missing subsurface
	  BSDF += albedo * blender_calculate_disney_brdf_diffuse_mul(si, roughness, bs, active);
	}
	BSDF = BSDF * diffuse_weight;
	//missing sheen

	if(any_or<true>(specular_weight > 1e-5f)){
	  //assume no anisotropic
	  Float aspect = 1.0f;
	  Float r2 = roughness * roughness;
	  Float alpha_x = r2 / aspect;
	  Float alpha_y = r2 * aspect;
	  //assume speculartint = 0.0
	  Spectrum tmp_col = Spectrum(1.0f) * (1.0f - specular_tint) + ctint * specular_tint;
	  //Spectrum Cspec0 = (specular * 0.08f * tmp_col) * (1.0f - metallic) + albedo * metallic; //blender implementation
	  //lerp(a,b,t) b * t + (1-t)*a or (1-t)*a+t*b
	  Spectrum Cspec0 = lerp(metallic, pbrt_schlickR0FromEta(ior) * tmp_col, albedo);//FIXME: we should use ior based implementation
	  MicrofacetDistribution distr(MicrofacetType::GGX, roughness, true);
	  Normal3f m;
	  Float pdf;
	  Throw("we need to reuse sample 2 here");
	  //WARNING: we need to reuse sample_2!!!!!! This is super bad since albedo direction we have two outgoing direction now. we need sample twice.
	  std::tie(m, pdf) = distr.sample(si.wi, sample2);
	  bs.wo = reflect(si.wi, m);
	  //actual disney fresnel uses mixture of dielectric  and schlick
	  Spectrum F =  blender_fresnel(Cspec0, ior, m, bs.wo);
	  Spectrum specular_comp = distr.eval(m) * distr.G(si.wi, bs.wo, m) * F; // missing fresnel term
	}
	
	// if (unlikely(none_or<false>(active) ||
	// 			 !ctx.is_enabled(BSDFFlags::DiffuseReflection)))
	//   return { bs, 0.f };
	
	
		
	bs.wo = warp::square_to_cosine_hemisphere(sample2);
	bs.pdf = warp::square_to_cosine_hemisphere_pdf(bs.wo);
	bs.eta = 1.f;
	bs.sampled_type = +BSDFFlags::DiffuseReflection;
	bs.sampled_component = 0;

	
	//select(active && bs.pdf > 0.f, value, 0.f)
	return { bs, BSDF };
  }
  
  Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
				const Vector3f &wo, Mask active) const override {
	// MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

	// if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
	//   return 0.f;

	// Float cos_theta_i = Frame3f::cos_theta(si.wi),
	//   cos_theta_o = Frame3f::cos_theta(wo);

	// active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

	// UnpolarizedSpectrum value =
	//   m_reflectance->eval(si, active) * math::InvPi<Float> * cos_theta_o;

	// return select(active, unpolarized<Spectrum>(value), 0.f);
	return 0.0f;
  }

  Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
			const Vector3f &wo, Mask active) const override {
	// MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

	// if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
	//   return 0.f;

	// Float cos_theta_i = Frame3f::cos_theta(si.wi),
	//   cos_theta_o = Frame3f::cos_theta(wo);

	// Float pdf = warp::square_to_cosine_hemisphere_pdf(wo);

	// return select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
	return 0.0f;
  }

  void traverse(TraversalCallback *callback) override {
	callback->put_object("reflectance", m_reflectance.get());
  }

  std::string to_string() const override {
	std::ostringstream oss;
	oss << "DisneyBRDF[" << std::endl
		<< "  reflectance = " << string::indent(m_reflectance) << std::endl
		<< "]";
	return oss.str();
  }

  MTS_DECLARE_CLASS()
  private:
  ref<Texture> m_reflectance;
  ref<Texture> m_metallic;
  ref<Texture> m_roughness;
  Float m_ior;
};

MTS_IMPLEMENT_CLASS_VARIANT(DisneyBRDF, BSDF)
MTS_EXPORT_PLUGIN(DisneyBRDF, "Disney BRDF Implementation")
NAMESPACE_END(mitsuba)
