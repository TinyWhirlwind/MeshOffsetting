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
		SlicerSDF(const CMeshO& mesh, Scalarm offset, int nodeCount = 129, int blockWidth = 8);
		~SlicerSDF() = default;
		SlicerSDF(const SlicerSDF&) = delete;
		SlicerSDF& operator=(const SlicerSDF&) = delete;

		void StartSlicerOne(int layer, Scalarm layer_thickness);
		void StartSlicerAll(int layer_count, Scalarm layer_thickness);

		SliceResult SliceAtZ(Scalarm zi);
		const std::vector<SliceResult>& GetLayers() const;


		//测试
		const std::vector<Point3m>& GetErrorNode();
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

		void BuildPlaneField(Scalarm zi);
		bool HasStoredNodeInCell(
			int x,
			int y,
			int lower_z,
			bool needs_upper);
		bool FillPlaneNode(
			int x,
			int y,
			int lower_z,
			Scalarm weight,
			bool needs_upper);
		void ProcessCell(int x, int y, SliceResult& result);
		std::size_t GetOrCreateCrossing(int ax, int ay, int bx, int by, Scalarm zi);
		void Connect(std::size_t a, std::size_t b);
		void TraceContours(SliceResult& result);
		static bool SamePoint(const Point3m& a, const Point3m& b);
		static bool HasArea(const std::vector<Point3m>& contour);
	private:
		NarrowBandSDF _nsdf;
		UniformGrid* _grid = nullptr;

		Scalarm _offset = 0.0;
		Scalarm _cell_width = 0;

		PlaneField _plane_grid;
		Box3m _bbox;

		int _size_x = 0;
		int _size_y = 0;
		int _size_z = 0;
		std::vector<ContourNode> _contour_nodes;

		// 网格边两端的顶点 ID -> 轮廓交点 ID。
		std::map<std::pair<std::size_t, std::size_t>, std::size_t> _crossing_ids;
		std::vector<SliceResult> _layers;
		std::vector<Point3m> _error_nodes;
	};
}
