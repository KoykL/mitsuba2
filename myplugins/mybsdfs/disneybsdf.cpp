#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/string.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)


template <typename Float, typename Spectrum>
class DisneyBSDF final : public BSDF<Float, Spectrum> {
public:
  MTS_IMPORT_BASE(BSDF, m_flags, m_components)
  MTS_IMPORT_TYPES(Texture)
protected:
  ref<Base> diffuse;
  ref<Base> specular;
  ref<Texture> m_metallic;
  
  ref<Base> get_bsdf_with_name(const std::string &name, const Properties &props) const{
	PluginManager *mgr = PluginManager::instance();
	const Class *plugin_class = mgr->get_plugin_class(name, MTS_CLASS(Base)->variant());
	
    Assert(plugin_class != nullptr);
    ref<Object> object = plugin_class->construct(props);
	Assert(object->class_()->derives_from(MTS_CLASS(Base)));
	return static_cast<Base *>(object.get());
  } 
public:
  
  DisneyBSDF(const Properties &props) : Base(props) {
	diffuse = get_bsdf_with_name("disneydiffusebsdf", props);
	//diffuse = PluginManager::instance()->create_object<Base>(diffuse_props);
	//Properties specular_props = Properties(props);
	//specular_props.set_plugin_name("disneyspecularbsdf");
	//specular = PluginManager::instance()->create_object<Base>(specular_props);
	specular = get_bsdf_with_name("disneyspecularbsdf", props);
	m_metallic = props.texture<Texture>("m_metallic", 0.0f);
	m_components.clear();
	for (size_t j = 0; j < diffuse->component_count(); ++j){
	  m_components.push_back(diffuse->flags(j));
	}
	for (size_t j = 0; j < diffuse->component_count(); ++j){
	  m_components.push_back(specular->flags(j));
	}
	
	m_flags = diffuse->flags() | specular->flags();
  }

  std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
										   const SurfaceInteraction3f &si,
										   Float sample1,
										   const Point2f &sample2,
										   Mask active) const override {
	MTS_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
	Float weight = eval_weight(si, active);
	if (unlikely(ctx.component != (uint32_t) -1)) {
	  bool sample_first = ctx.component < diffuse->component_count();
	  BSDFContext ctx2(ctx);
	  if (!sample_first){
		ctx2.component -= (uint32_t) diffuse->component_count();
	  }
	  else{
		weight = 1.f - weight;
	  }
	  std::pair<BSDFSample3f, Spectrum> result;
	  if (sample_first){
		result = diffuse->sample(ctx2, si, sample1, sample2, active);
	  }else{
		result = specular->sample(ctx2, si, sample1, sample2, active);
	  }
	  
	  result.second *= weight;
	  return result;
	}

	BSDFSample3f bs = zero<BSDFSample3f>();
	Spectrum result(0.f);

	Mask m0 = active && sample1 >  weight,
	  m1 = active && sample1 <= weight;

	if (any_or<true>(m0)) {
	  auto [bs0, result0] = diffuse->sample(ctx, si, (sample1 - weight) / (1 - weight), sample2, m0);
	  masked(bs, m0) = bs0;
	  masked(result, m0) = result0;
	}

	if (any_or<true>(m1)) {
	  auto [bs1, result1] = specular->sample(ctx, si, sample1 / weight, sample2, m1);
	  masked(bs, m1) = bs1;
	  masked(result, m1) = result1;
	}

	return { bs, result };
  }

  Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
				const Vector3f &wo, Mask active) const override {
	MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

	Float weight = eval_weight(si, active);
	if (unlikely(ctx.component != (uint32_t) -1)) {
	  bool sample_first = ctx.component < diffuse->component_count();
	  BSDFContext ctx2(ctx);
	  if (!sample_first){
		ctx2.component -= (uint32_t) diffuse->component_count();
	  }
	  else{
		weight = 1.f - weight;
	  }
	  Spectrum result;
	  if (sample_first){
		result = diffuse->eval(ctx2, si, wo, active);
	  }else{
		result = specular->eval(ctx2, si, wo, active);
	  }
	  return weight * result;
	}
	
	return diffuse->eval(ctx, si, wo, active) * (1 - weight) + specular->eval(ctx, si, wo, active) * weight;
  }

  Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
			const Vector3f &wo, Mask active) const override {
	MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

	if (unlikely(ctx.component != (uint32_t) -1)) {
	  bool sample_first = ctx.component < diffuse->component_count();
	  BSDFContext ctx2(ctx);
	  if (!sample_first){
		ctx2.component -= (uint32_t) diffuse->component_count();
	  }
	  Float pdf;
	  if(sample_first){
		pdf = diffuse->pdf(ctx2, si, wo, active);
	  }else{
		pdf = specular->pdf(ctx2, si, wo, active);
	  }
	  return pdf;
	}
	
	Float weight = eval_weight(si, active);
	return diffuse->pdf(ctx, si, wo, active) * (1 - weight) + specular->pdf(ctx, si, wo, active) * weight;
  }

  Float eval_weight(const SurfaceInteraction3f &si, Mask active) const{
	Float metallic = m_metallic->eval_1(si, active);
	Float diffuse_weight = (1.0f - clamp(metallic, 0.0f, 1.0f)); //* (1.0 - calmp(transmissionm 0.0, 1.0))
	Float specular_weight = 1.0;// - final_transmission
	Float weight = specular_weight / (diffuse_weight + specular_weight);
	return weight;
  }
  
  void traverse(TraversalCallback *callback) override {
	callback->put_object("diffuse", diffuse.get());
	callback->put_object("specular", specular.get());
  }
  
  std::string to_string() const override {
	std::ostringstream oss;
	oss << "DisneyBSDF[" << std::endl
		<< "  diffuse = " << string::indent(diffuse) << "," << std::endl
		<< "  specular = " << string::indent(specular) << std::endl
		<< "]";
	return oss.str();
  }
  
  MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_VARIANT(DisneyBSDF, BSDF)
MTS_EXPORT_PLUGIN(DisneyBSDF, "DisneyBSDF material")
NAMESPACE_END(mitsuba)
