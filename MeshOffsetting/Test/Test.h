/**
 * @file TestSlicer.h
 * @author zhangpingping
 * @date 2026-11-19
 * @brief
 *
 */
#ifndef MESHOFFSETTING_TEST_H
#define MESHOFFSETTING_TEST_H

#include "base_type.h"

#include <utility>
#include <vector>
namespace meshoffsetting
{
	class TestSlicer
	{
	public:
		TestSlicer();
		~TestSlicer() = default;
		TestSlicer(const TestSlicer&) = delete;
		TestSlicer& operator=(const TestSlicer&) = delete;

		void SetSize(int x,int y,int z);
		void GetLayerNode(int layer);
		void SetLayerNode(int layer, const std::vector<Point3m>& nodes);
	private:
		int _size_x = 0;
		int _size_y = 0;
		int _size_z = 0;
		std::vector<std::pair<int, std::vector<Point3m>>> xoy_node_list;
	};
}

#endif
