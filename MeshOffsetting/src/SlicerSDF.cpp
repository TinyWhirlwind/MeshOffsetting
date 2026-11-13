#include "SlicerSDF.h"
#include "Tolerance.h"
namespace meshoffsetting
{
	SlicerSDF::SlicerSDF(const CMeshO& mesh, Scalarm offset):
		_nsdf(mesh, offset),
		_offset(offset)
	{
		_grid = _nsdf.GetGird();
		_bbox = _grid->GetBox();
		_plane_grid.resize(_grid->GetSize().X(), _grid->GetSize().Y());
	}

	void SlicerSDF::StartSlicerOne(int layer, Scalarm layer_thickness)
	{
		Scalarm z = (layer + 1) * layer_thickness;

	}

	void SlicerSDF::StartSlicerAll(int layer_count, Scalarm layer_thickness)
	{
		Point3i grid_size = _grid->GetSize();
		Scalarm cell_width = _grid->GetCellWidth();
		Scalarm z_min = _bbox.min.Z();
		Scalarm z_max = _bbox.max.Z();
		for (int layer = 0; layer < layer_count; ++layer)
		{
			Scalarm zi = (layer + 1) * layer_thickness;
			if (zi < z_min || zi > z_max) 
			{
				return;
			}

			_plane_grid.clear();
			_segments.clear();
			Scalarm real_z = (zi - z_min) / cell_width;
			int integer_z = static_cast<int>(std::floor(real_z));
			Scalarm extra_z = real_z - integer_z;
			if (extra_z * cell_width > kScalarTolerance)
			{
				//位于两层之间
				for (int x = 0; x < grid_size.X(); ++x)
				{
					for (int y = 0; y < grid_size.Y(); ++y)
					{
						GridNode& nd0 = _grid->Node(Point3i(x, y, integer_z));
						GridNode& nd1 = _grid->Node(Point3i(x, y, integer_z + 1));
						if (nd0.state != GridNodeState::NodeDistanceReady
							|| nd1.state != GridNodeState::NodeDistanceReady)
						{
							continue;
						}
						Scalarm slicer_sdf = nd0.signed_distance * (1 - extra_z) + nd1.signed_distance* extra_z;
						int node_id = x * grid_size.X() + y;
						_plane_grid.values[x][y] = slicer_sdf - _offset;
						_plane_grid.valid[x][y] = true;
					}
				}
			}
			else
			{
				//正好经过一层 Grid Node
				for (int x = 0; x < grid_size.X(); ++x)
				{
					for (int y = 0; y < grid_size.Y(); ++y)
					{
						GridNode& nd = _grid->Node(Point3i(x, y, integer_z));
						if (nd.state == GridNodeState::NodeInvalid)
						{
							continue;
						}
						_plane_grid.values[x][y] = nd.signed_distance - _offset;
						_plane_grid.valid[x][y] = true;
					}
				}
			}

		}
	}

	void SlicerSDF::ApplySlicer(Scalarm slicer_z)
	{
		for (int i = 0; i < _plane_grid.values.size() - 1; ++i)
		{
			for (int j = 0; j < _plane_grid.values[i].size() - 1; ++j)
			{
				GetSgementList(i, j, slicer_z);
			}
		}

		//查找闭合轮廓

	}

	int SlicerSDF::GetCaseId(bool n00, bool n01, bool n11, bool n10)
	{
		if ( n00 &&  n10 &&  n11 &&  n01)return 15;
		if ( n00 && !n10 &&  n11 &&  n01)return 14;
		if ( n00 &&  n10 && !n11 &&  n01)return 13;
		if ( n00 && !n10 && !n11 &&  n01)return 12;
		if ( n00 &&  n10 &&  n11 && !n01)return 11;
		if ( n00 && !n10 &&  n11 && !n01)return 10;
		if ( n00 &&  n10 && !n11 && !n01)return 9;
		if ( n00 && !n10 && !n11 && !n01)return 8;
		if (!n00 &&  n10 &&  n11 &&  n01)return 7;
		if (!n00 && !n10 &&  n11 &&  n01)return 6;
		if (!n00 &&  n10 && !n11 &&  n01)return 5;
		if (!n00 && !n10 && !n11 &&  n01)return 4;
		if (!n00 &&  n10 &&  n11 && !n01)return 3;
		if (!n00 && !n10 &&  n11 && !n01)return 2;
		if (!n00 &&  n10 && !n11 && !n01)return 1;
		if (!n00 && !n10 && !n11 && !n01)return 0;
		return -1;
	}

	void SlicerSDF::GetSgementList(int i,int j, Scalarm slicer_z)
	{
		bool valid[4] = {};
		Scalarm values[4] = {};
		Point3m node_postions[4] = {};
		std::pair<int, int> node_pairs[4] = {};
		Scalarm cell_width = _grid->GetCellWidth();

		int order = 0;
		for (int m = 0; m <= 1; ++m)
		{
			for (int n = 0; n <= 1; ++n)
			{
				if (order == 2)
				{
					order = 3;
					valid[order] = _plane_grid.valid[i + m][j + n];
					values[order] = _plane_grid.values[i + m][j + n];
					node_postions[order] = Point3m{ _bbox.min.X() + i * cell_width,
												_bbox.min.Y() + j * cell_width,
												slicer_z };	
					node_pairs[order] = std::make_pair(i, j);
					--order;
				}
				else
				{
					++order;
				}
			}
		}
		/*
		v0**vi0****v1
		 **********
		vi3*******vi1
		 **********
		v3***vi2***v2
		*/
		auto ComputerInterpolation = [values, node_postions, node_pairs](int edge_index) ->ContourPoint
			{
				Scalarm v0 = values[edge_index];
				Scalarm v1 = values[(edge_index + 1) % 4];
				bool crossed = (v0 < 0) != (v1 < 0);
				Scalarm t = v0 / (v0 - v1);
				Point3m p0 = node_postions[edge_index];
				Point3m p1 = node_postions[(edge_index + 1) % 4];
				Point3m interpolation_p = p0 + (p1 - p0) * t;
				return ContourPoint(node_pairs->first, node_pairs->second, interpolation_p);
			};

		std::uint8_t case_id = GetCaseId(valid[0], valid[1], valid[2], valid[3]);
		switch (case_id)
		{
		case 0:
		case 15:
			break;
		case 5:
		case 10:
			//对角线同号，属于二义性情况
			break;
		case 1:
		{
			ContourPoint vi2 = ComputerInterpolation(2);
			ContourPoint vi3 = ComputerInterpolation(3);
			_segments.push_back(ContourSegment(vi2, vi3));
			break;
		}
			/*
		case 2:
			break;
		case 3:
			break;
		case 4:
			break;*/
		default:
			break;
		}
	}
}