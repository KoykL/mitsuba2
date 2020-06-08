#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/microfacet.h>
#include <mitsuba/render/fresnel.h>
#include <iostream>
NAMESPACE_BEGIN(mitsuba)
template <typename Float, typename Spectrum>
class DisneyDiffuseBRDF final : public BSDF<Float, Spectrum> {
public:
  MTS_IMPORT_BASE(BSDF, m_flags, m_components)
  MTS_IMPORT_TYPES(Texture, MicrofacetDistribution)

  DisneyDiffuseBRDF(const Properties &props) : Base(props) {
	m_reflectance = props.texture<Texture>("reflectance", .5f);
	m_roughness = props.texture<Texture>("roughness", 1.0f);
	m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
	m_components.push_back(m_flags);
  }
  Float blender_schlick_fresnel(Float u) const{
	Float m = clamp(1.0f - u, 0.0f, 1.0f);
	Float m2 = m * m;
	return m2 * m2 * m;
  }
  Float blender_calculate_disney_brdf_diffuse_mul(const SurfaceInteraction3f &si, const Float &roughness, const Vector3f &wo, Mask &active) const{
	//Vector3f world_L = si.sh_frame.to_world(wo);
	//Float NgdotL = dot(si.n, world_L);
	//active &= NgdotL > 0.f;
	//std::cout << si.n << std::endl;
  	Float NdotL = Frame3f::cos_theta(wo); //outcome
	active &= NdotL > 0.f;
	Float NdotV = Frame3f::cos_theta(si.wi); //income
	active &= NdotV > 0.f;
	if (unlikely(none_or<false>(active))){
	  return 0.f;
	}
	
	Vector3f H = normalize(si.wi + wo);
	Float LdotH = dot(wo, H);
	Float FL = blender_schlick_fresnel(NdotL);
	Float FV = blender_schlick_fresnel(NdotV);
	Float Fd90 = 0.5f + 2.0f * LdotH * LdotH * roughness;
	//lerp(a,b,t) b * t + (1-t)*a
	Float Fd = (1.0f * (1.0f - FL) + Fd90 * FL) * (1.0f * (1.0f - FV) + Fd90 * FV);
	Float value = 1.0f / Float(M_PI) * NdotL * Fd; //cos here
	//std::cout << "still active" << active << std::endl;
	return select(active, value, 0.0f);
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
	BSDFSample3f bs = zero<BSDFSample3f>();
	if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
	  return {bs, 0.f};

	//We might need to use packet!!!
	Spectrum albedo = m_reflectance->eval(si, active);
	Float roughness = m_roughness->eval_1(si, active);
	
	Spectrum BSDF(0.0f);
		
	bs.wo = warp::square_to_cosine_hemisphere(sample2);
	bs.pdf = warp::square_to_cosine_hemisphere_pdf(bs.wo);
	
	BSDF = albedo * blender_calculate_disney_brdf_diffuse_mul(si, roughness, bs.wo, active);
	//cout << "reached_here" << endl;
	// if (unlikely(none_or<false>(active) ||
	// 			 !ctx.is_enabled(BSDFFlags::DiffuseReflection)))
	//   return { bs, 0.f };		

	bs.eta = 1.f;
	bs.sampled_type = +BSDFFlags::DiffuseReflection;
	bs.sampled_component = 0;
	//std::cout << BSDF << std::endl;
	//cos is multiplied in blender_calculate_disney_brdf_diffuse_mul
	return { bs, select(active && bs.pdf > 0.f, BSDF/bs.pdf, 0.f) };
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
	MTS_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
	if constexpr (is_polarized_v<Spectrum>) {
        Throw("DisneyBRDF: I don't know how to support polarized mode");
    }
	BSDFSample3f bs = zero<BSDFSample3f>();
	if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
	  return 0.f;

	Spectrum albedo = m_reflectance->eval(si, active);
	Float roughness = m_roughness->eval_1(si, active);
	
	Spectrum BSDF = albedo * blender_calculate_disney_brdf_diffuse_mul(si, roughness, wo, active);
	//cout << "reached_here" << endl;
	// if (unlikely(none_or<false>(active) ||
	// 			 !ctx.is_enabled(BSDFFlags::DiffuseReflection)))
	//   return { bs, 0.f };		
	return BSDF;
  }

  Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
			const Vector3f &wo, Mask active) const override {
	MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

	// if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
	//   return 0.f;

	// Float cos_theta_i = Frame3f::cos_theta(si.wi),
	//   cos_theta_o = Frame3f::cos_theta(wo);

	// Float pdf = warp::square_to_cosine_hemisphere_pdf(wo);

	// return select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
	if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
	  return 0.f;
	// Vector3f world_L = si.sh_frame.to_world(wo);
	// Float NgdotL = dot(si.n, world_L);
	// active &= NgdotL > 0.f;
  	Float NdotL = Frame3f::cos_theta(wo); //outcome
	active &= NdotL > 0.f;
	Float NdotV = Frame3f::cos_theta(si.wi); //income
	active &= NdotV > 0.f;

	Float pdf = warp::square_to_cosine_hemisphere_pdf(wo);
	return select(active, pdf, 0.0f);
  }

  void traverse(TraversalCallback *callback) override {
	callback->put_object("reflectance", m_reflectance.get());
  }

  std::string to_string() const override {
	std::ostringstream oss;
	oss << "DisneyDiffuseBRDF[" << std::endl
		<< "  reflectance = " << string::indent(m_reflectance) << std::endl
		<< "]";
	return oss.str();
  }

  MTS_DECLARE_CLASS()
  private:
  ref<Texture> m_reflectance;
  ref<Texture> m_roughness;
  Float m_ior;
};

MTS_IMPLEMENT_CLASS_VARIANT(DisneyDiffuseBRDF, BSDF)
MTS_EXPORT_PLUGIN(DisneyDiffuseBRDF, "Disney Diffuse Implementation")
NAMESPACE_END(mitsuba)
