/**
 * @file SlicerSDF.h
 * @author zhangpingping
 * @date 2026-11-12
 * @brief
 *
 */
#include "NarrowBandSDF.h"
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
			std::vector<std::vector<Scalarm>> values;
			std::vector<std::vector<bool>> valid;
		};

		// 轮廓交点，不是原来的网格采样节点。
		struct ContourNode
		{
			Point3m point;
			std::vector<std::size_t> neighbours;
		};

		struct LayerWorkspace
		{
			PlaneField plane_grid;
			std::vector<ContourNode> contour_nodes;
			// 网格边两端的顶点 ID -> 轮廓交点 ID。
			std::map<std::pair<std::size_t, std::size_t>, std::size_t> crossing_ids;
			std::vector<Point3m> error_nodes;
		};

		void InitializeWorkspace(LayerWorkspace& workspace) const;
		bool ComputeSlicePosition(
			Scalarm zi,
			int& lower_z,
			Scalarm& weight,
			bool& needs_upper) const;
		bool PrepareDistanceCache(const std::vector<Scalarm>& slice_z_values);
		bool GetGridNodeDistance(
			int x,
			int y,
			int z,
			Scalarm& distance);
		bool GetCachedGridNodeDistance(
			int x,
			int y,
			int z,
			Scalarm& distance) const;
		SliceResult SliceAtZImpl(Scalarm zi, LayerWorkspace& workspace) const;
		bool BuildPlaneField(Scalarm zi, LayerWorkspace& workspace) const;
		bool FillPlaneNode(
			int x,
			int y,
			int lower_z,
			Scalarm weight,
			bool needs_upper,
			LayerWorkspace& workspace) const;
		void ProcessCell(int x, int y, SliceResult& result, LayerWorkspace& workspace) const;
		std::size_t GetOrCreateCrossing(
			int ax,
			int ay,
			int bx,
			int by,
			Scalarm zi,
			LayerWorkspace& workspace) const;
		void Connect(std::size_t a, std::size_t b, LayerWorkspace& workspace) const;
		void TraceContours(SliceResult& result, LayerWorkspace& workspace) const;
		static bool SamePoint(const Point3m& a, const Point3m& b);
		static bool HasArea(const std::vector<Point3m>& contour);
	private:
		NarrowBandSDF _nsdf;
		UniformGrid* _grid = nullptr;

		Scalarm _offset = 0.0;
		Scalarm _cell_width = 0;

		std::vector<Scalarm> _distance_cache;
		std::vector<unsigned char> _distance_cached;
		Box3m _bbox;

		int _size_x = 0;
		int _size_y = 0;
		int _size_z = 0;
		std::vector<SliceResult> _layers;
		//std::vector<Point3m> _error_nodes;
	};
}
