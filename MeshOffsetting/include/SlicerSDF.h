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
	class SlicerSDF
	{
		struct PlaneField
		{
			std::vector<std::vector<Scalarm>> values;
			std::vector<std::vector<bool>> valid;
			void clear()
			{
				for (auto& row : values)
				{
					std::fill(row.begin(), row.end(), Scalarm(0));
				}

				for (auto& row : valid)
				{
					std::fill(row.begin(), row.end(), false);
				}
			}

			void resize(int size_x, int size_y)
			{
				values.assign(
					size_x,
					std::vector<Scalarm>(size_y, Scalarm(0)));
				valid.assign(
					size_x,
					std::vector<bool>(size_y, false));
			}
		};

		struct ContourPoint
		{
			int begin = -1;
			int end = -1;
			Point3m point{0, 0, 0};
			ContourPoint(int v0, int v1, Point3m p)
			{
				begin = std::min(v0, v1);
				end = std::max(v0, v1);
				point = p;
			}
		};

		struct ContourSegment
		{
			ContourPoint start;
			ContourPoint end;

			ContourSegment(ContourPoint s, ContourPoint e):
				start(s), end(e)
			{
			}
		};
	public:
		SlicerSDF(const CMeshO& mesh, Scalarm offset);
		~SlicerSDF() = default;

		void StartSlicerOne(int layer, Scalarm layer_thickness);
		void StartSlicerAll(int layer_count, Scalarm layer_thickness);
	private:
		void ApplySlicer(Scalarm slicer_z);
		int MarchingSquares(bool n00, bool n10, bool n11, bool n01);
		int GetCaseId(bool n00, bool n10, bool n11, bool n01);
		void GetSgementList(int i, int j, Scalarm slicer_z);
	private:
		NarrowBandSDF _nsdf;
		Scalarm _offset = 0.0;
		UniformGrid* _grid = nullptr;

		PlaneField _plane_grid;
		Box3m _bbox;
		std::vector<ContourSegment> _segments;
	};
}