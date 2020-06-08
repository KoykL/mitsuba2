#include <vector>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/string.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#define TINYOBJLOADER_IMPLEMENTATION
#include <filesystem>
#include "tiny_obj_loader.h"
NAMESPACE_BEGIN(mitsuba)


template <typename Float, typename Spectrum>
class MultiMaterialBSDF final : public BSDF<Float, Spectrum> {
public:
  MTS_IMPORT_BASE(BSDF, m_flags, m_components)
  MTS_IMPORT_TYPES(Texture)
protected:
  std::vector<ref<Base>> m_bsdfs;
  std::vector<int> m_material_ids;
  ref<Base> get_bsdf_with_name(const std::string &name, const Properties &props) const{
	PluginManager *mgr = PluginManager::instance();
	const Class *plugin_class = mgr->get_plugin_class(name, MTS_CLASS(Base)->variant());
	
    Assert(plugin_class != nullptr);
    ref<Object> object = plugin_class->construct(props);
	Assert(object->class_()->derives_from(MTS_CLASS(Base)));
	return static_cast<Base *>(object.get());
  } 

public:
  MultiMaterialBSDF(const Properties &props) : Base(props) {
	int bsdf_index = 0;
	std::string path(props.string("filename"));
	
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	
	std::string warn;
	std::string err;
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, path.c_str());

	if (!warn.empty()) {
	  std::cout << warn << std::endl;
	}

	if (!err.empty()) {
	  std::cerr << err << std::endl;
	}

	if (!ret) {
	  Throw("MultiMaterialBSDF: tinyobjloader failed to parse the obj file");
	}
	m_material_ids.clear();
	for(size_t s = 0; s < shapes.size(); s++){
	  size_t index_offset = 0;
	  for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++){
		m_material_ids.push_back(shapes[s].mesh.material_ids[f]); // per face material ids
	  }
	}
	size_t num_materials = materials.size();
	for(size_t i = 0; i < num_materials; i++){
	  const ref<Object> & obj = props.object("bsdf_"+std::to_string(i));
	  const auto *cbsdf = dynamic_cast<const Base *>(obj.get());
	  auto *bsdf = const_cast<Base *>(cbsdf);
	  if(!bsdf){
		Throw("only bsdf can be specified");
	  }
	  m_bsdfs.push_back(bsdf);
	}
	m_components.clear();
	m_flags = 0;
	
	for (size_t i = 0; i < num_materials; ++i){
	  for (size_t j = 0; j < m_bsdfs[i]->component_count(); ++j){
		m_components.push_back(m_bsdfs[i]->flags(j));
	  }
	  m_flags |= m_bsdfs[i]->flags();
	}
  }
  
  
  std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
										   const SurfaceInteraction3f &si,
										   Float sample1,
										   const Point2f &sample2,
										   Mask active) const override {
	MTS_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
	typename SurfaceInteraction3f::Index face_idx = si.prim_index;	
	return m_bsdfs[face_idx]->sample(ctx, si, sample1, sample2, active);
  }

  Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
				const Vector3f &wo, Mask active) const override {
	MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);
	typename SurfaceInteraction3f::Index face_idx = si.prim_index;	
	return m_bsdfs[face_idx]->eval(ctx, si, wo, active);
  }

  Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
			const Vector3f &wo, Mask active) const override {
	MTS_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);
	typename SurfaceInteraction3f::Index face_idx = si.prim_index;	
	return m_bsdfs[face_idx]->pdf(ctx, si, wo, active);
  }
  
  void traverse(TraversalCallback *callback) override {
	for(size_t i = 0; i < m_bsdfs.size(); i++){
	  callback->put_object("bsdf_"+std::to_string(i), m_bsdfs[i].get());
	}
  }
  
  std::string to_string() const override {
	std::ostringstream oss;
	oss << "MultiMaterialBSDF[" << std::endl
		<< "]";
	return oss.str();
  }
  
  MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_VARIANT(MultiMaterialBSDF, BSDF)
MTS_EXPORT_PLUGIN(MultiMaterialBSDF, "MultiMaterialBSDF material")
NAMESPACE_END(mitsuba)
