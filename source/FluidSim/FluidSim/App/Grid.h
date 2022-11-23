#pragma once

#include <vector>
#include <memory>

#include "Domain.h"
#include "GridNode.h"

class Grid
{
public:
	Grid() = default;
	~Grid() = default;

	Grid(float inLength, float inWidth, float inHeight, float inCellSize);
	Grid(const Domain& inDomain, float inCellSize);

	void StepGrid(float deltaTime);

	GridNode& GetGridNode(int index) { return mGridNodes[index]; }
	std::vector<GridNode>& GetGridNodes() { return mGridNodes; }
	int GetNumGridNodes() { return mGridNodes.size(); }

private:

	void InitializeGridNodes();

	float mLength;
	float mWidth;
	float mHeight;

	float mCellSize;

	// 'Double buffer' grid nodes
	std::vector<GridNode> mGridNodes;
	std::vector<GridNode> mGridNodesTemp;
};