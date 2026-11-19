/**
 * @file SlicerSDF.h
 * @author zhangpingping
 * @date 2026-11-12
 * @brief
 *
 */
#ifndef MESHOFFSETTING_SLICER_SDF_H
#define MESHOFFSETTING_SLICER_SDF_H

#include "NarrowBandSDF.h"
#include "Test/Test.h"
namespace meshoffsetting
{
	struct SliceResult
	{
		Scalarm z = 0;
		std::vector<std::vector<Point3m>> contours;
		std::size_t skipped_cells = 0;
	};

	class SlicerSDF
	{
	public:
		SlicerSDF(const CMeshO& mesh, Scalarm offset, int nodeCount = 257, int blockWidth = 8);
		~SlicerSDF() = default;
		SlicerSDF(const SlicerSDF&) = delete;
		SlicerSDF& operator=(const SlicerSDF&) = delete;

		void StartSlicerOne(int layer, Scalarm layer_thickness);
		void StartSlicerAll(int layer_count, Scalarm layer_thickness);

		SliceResult SliceAtZ(Scalarm zi);
		const std::vector<SliceResult>& GetLayers() const;


		//测试
		//const std::vector<Point3m>& GetErrorNode();
	private:
		struct PlaneField
		{
			// Flat node caches; only touched nodes are reset between slices.
			std::vector<Scalarm> values;
			std::vector<unsigned char> valid;
			// Reuse candidates while the Z planes and interpolation mode match.
			std::vector<std::size_t> candidate_cell_indices;
			std::vector<std::size_t> touched_nodes;
			int candidate_lower_z = -1;
			bool candidate_needs_upper = false;
		};

		enum class PlaneNodeState : unsigned char
		{
			Missing,
			NarrowBand,
			Support,
			Failed
		};

		struct DistancePlane
		{
			int z = -1;
			std::vector<Scalarm> values;
			std::vector<PlaneNodeState> states;
		};

		// 轮廓交点，不是原来的网格采样节点。
		struct ContourNode
		{
			Point3m point;
			std::vector<std::size_t> neighbours;
		};

		struct LayerWorkspace
		{
			DistancePlane lower_plane;
			DistancePlane upper_plane;
			PlaneField plane_grid;
			std::vector<ContourNode> contour_nodes;
			// 网格边两端的顶点 ID -> 轮廓交点 ID。
			std::map<std::pair<std::size_t, std::size_t>, std::size_t> crossing_ids;
			std::vector<Point3m> error_nodes;
		};

		void InitializeWorkspace(LayerWorkspace& workspace) const;
		bool ComputeSlicePosition(Scalarm zi, int& lower_z, Scalarm& weight, bool& needs_upper) const;
		bool LoadDistancePlane(int z, DistancePlane& plane) const;
		bool PrepareDistancePlanes(int lower_z,bool needs_upper,LayerWorkspace& workspace) const;
		bool IsCandidateCell(int x,int y,bool needs_upper,const LayerWorkspace& workspace) const;
		bool EnsurePlaneNodeDistance(int x,int y,DistancePlane& plane) const;
		std::size_t PlaneIndex(int x, int y) const;
		SliceResult SliceAtZImpl(Scalarm zi, LayerWorkspace& workspace) const;
		bool BuildPlaneField(Scalarm zi, LayerWorkspace& workspace) const;
		bool FillPlaneNode(int x,int y,Scalarm weight,bool needs_upper,LayerWorkspace& workspace) const;
		void ProcessCell(int x, int y, SliceResult& result, LayerWorkspace& workspace) const;
		std::size_t GetOrCreateCrossing(int ax,int ay,int bx,int by,Scalarm zi,LayerWorkspace& workspace) const;
		void Connect(std::size_t a, std::size_t b, LayerWorkspace& workspace) const;
		void TraceContours(SliceResult& result, LayerWorkspace& workspace) const;
		static bool SamePoint(const Point3m& a, const Point3m& b);
		static bool HasArea(const std::vector<Point3m>& contour);
	private:
		NarrowBandSDF _nsdf;
		TestSlicer _test;

		Scalarm _offset = 0.0;
		Scalarm _cell_width = 0;

		Box3m _bbox;

		int _size_x = 0;
		int _size_y = 0;
		int _size_z = 0;
		std::vector<SliceResult> _layers;
		//std::vector<Point3m> _error_nodes;
	};
}

#endif
