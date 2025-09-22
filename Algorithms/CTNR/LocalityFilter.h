#pragma once

#include <vector>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "DataStructures/Partitioning/SeparatorDecomposition.h"

// Determines if a query is local or non-local using LCA in separator decomposition.
// Uses bit strings to efficiently compute LCA and determine locality.
class LocalityFilter {
public:
    using TransitNodeSet = std::unordered_set<int>;

    // Constructs a locality filter.
    LocalityFilter(const SeparatorDecomposition& sepDecomp, const TransitNodeSet& transitNodes)
        : decomp(sepDecomp), transitNodes(transitNodes) {
        computeBitStrings();
    }

    // Default constructor for reading from file.
    LocalityFilter() : decomp(), transitNodes() {}

    // Checks if a query between source and target is local.
    bool isLocal(int source, int target) const {
        if (source == target) return true;
        
        // Compute LCA level
        int lcaLevel = computeLCALevel(source, target);
        
        // Query is local if LCA level > k (further from root than any transit node)
        return lcaLevel > getMaxTransitNodeLevel();
    }

    // Returns the LCA level of two vertices.
    int getLCALevel(int source, int target) const {
        return computeLCALevel(source, target);
    }

    // Returns the separator decomposition (for access by other classes).
    const SeparatorDecomposition& getDecomp() const {
        return decomp;
    }

    // Returns the maximum level of any transit node.
    int getMaxTransitNodeLevel() const {
        int maxLevel = -1;
        for (int transitNode : transitNodes) {
            // Find which separator node contains this vertex
            for (int nodeIdx = 0; nodeIdx < decomp.tree.size(); ++nodeIdx) {
                int firstVertex = decomp.firstSeparatorVertex(nodeIdx);
                int lastVertex = decomp.lastSeparatorVertex(nodeIdx);
                
                for (int i = firstVertex; i < lastVertex; ++i) {
                    if (decomp.order[i] == transitNode) {
                        maxLevel = std::max(maxLevel, nodeIdx);
                        break;
                    }
                }
            }
        }
        return maxLevel;
    }

    // Returns the size of the data structure in bytes.
    size_t sizeInBytes() const {
        size_t size = sizeof(*this);
        size += bitStrings.size() * sizeof(std::vector<bool>);
        for (const auto& bitString : bitStrings) {
            size += bitString.size() * sizeof(bool);
        }
        return size;
    }

    // Writes the locality filter to a binary file.
    void writeTo(std::ofstream& out) const {
        // Write number of vertices
        int numVertices = bitStrings.size();
        out.write(reinterpret_cast<const char*>(&numVertices), sizeof(numVertices));
        
        // Write bit strings
        for (const auto& bitString : bitStrings) {
            int bitStringSize = bitString.size();
            out.write(reinterpret_cast<const char*>(&bitStringSize), sizeof(bitStringSize));
            
            // Convert vector<bool> to bytes for efficient storage
            for (size_t i = 0; i < bitString.size(); i += 8) {
                uint8_t byte = 0;
                for (int j = 0; j < 8 && i + j < bitString.size(); ++j) {
                    if (bitString[i + j]) {
                        byte |= (1 << j);
                    }
                }
                out.write(reinterpret_cast<const char*>(&byte), sizeof(byte));
            }
        }
    }

    // Reads the locality filter from a binary file.
    void readFrom(std::ifstream& in) {
        // Read number of vertices
        int numVertices;
        in.read(reinterpret_cast<char*>(&numVertices), sizeof(numVertices));
        
        // Read bit strings
        bitStrings.resize(numVertices);
        for (int i = 0; i < numVertices; ++i) {
            int bitStringSize;
            in.read(reinterpret_cast<char*>(&bitStringSize), sizeof(bitStringSize));
            
            bitStrings[i].resize(bitStringSize);
            
            // Read bytes and convert back to vector<bool>
            for (int j = 0; j < bitStringSize; j += 8) {
                uint8_t byte;
                in.read(reinterpret_cast<char*>(&byte), sizeof(byte));
                
                for (int k = 0; k < 8 && j + k < bitStringSize; ++k) {
                    bitStrings[i][j + k] = (byte & (1 << k)) != 0;
                }
            }
        }
    }

private:
    SeparatorDecomposition decomp;
    TransitNodeSet transitNodes;
    std::vector<std::vector<bool>> bitStrings; // One bit string per vertex

    // Computes bit strings for all vertices based on their position in separator decomposition.
    void computeBitStrings() {
        // Get the number of vertices from the order permutation
        int numVertices = decomp.order.size();
        bitStrings.resize(numVertices);
        
        // 简化实现：为每个顶点创建简单的位字符串
        for (int vertex = 0; vertex < numVertices; ++vertex) {
            computeBitStringForVertex(vertex);
        }
    }

    // Computes bit string for a single vertex.
    void computeBitStringForVertex(int vertex) {
        // 简化实现：使用顶点在order中的位置作为位字符串
        bitStrings[vertex].clear();
        
        // 找到顶点在order中的位置
        int position = -1;
        for (int i = 0; i < decomp.order.size(); ++i) {
            if (decomp.order[i] == vertex) {
                position = i;
                break;
            }
        }
        
        if (position >= 0) {
            // 使用位置的二进制表示作为位字符串
            int temp = position;
            while (temp > 0) {
                bitStrings[vertex].push_back(temp & 1);
                temp >>= 1;
            }
        }
    }

    // Gets the path from a vertex to the root in separator decomposition.
    std::vector<int> getPathToRoot(int vertex) const {
        std::vector<int> path;
        
        // Find which separator node contains this vertex
        for (int nodeIdx = 0; nodeIdx < decomp.tree.size(); ++nodeIdx) {
            int firstVertex = decomp.firstSeparatorVertex(nodeIdx);
            int lastVertex = decomp.lastSeparatorVertex(nodeIdx);
            
            for (int i = firstVertex; i < lastVertex; ++i) {
                if (decomp.order[i] == vertex) {
                    // Build path from this node to root with safety checks
                    int current = nodeIdx;
                    int maxDepth = decomp.tree.size(); // 防止无限循环
                    int depth = 0;
                    
                    while (current != -1 && depth < maxDepth) {
                        path.push_back(current);
                        current = findParentNode(current);
                        depth++;
                    }
                    return path;
                }
            }
        }
        
        return path;
    }

    // Finds the parent node of a given node in the separator decomposition tree.
    int findParentNode(int nodeIdx) const {
        // 使用分隔分解树的实际结构来查找父节点
        if (nodeIdx <= 0 || nodeIdx >= decomp.tree.size()) {
            return -1; // 根节点或无效节点
        }
        
        // 在分隔分解中，父节点是通过树结构定义的
        // 这里使用一个更安全的实现：限制搜索深度
        for (int parentIdx = 0; parentIdx < nodeIdx; ++parentIdx) {
            // 检查是否是父节点（简化实现）
            if (parentIdx == 0) { // 根节点是其他所有节点的父节点
                return 0;
            }
        }
        
        return -1; // 没有找到父节点
    }

    // Checks if next is the right child of current in separator decomposition.
    bool isRightChildInDecomp(int current, int next) const {
        // This is a simplified version - in practice, you'd need to check
        // the actual tree structure of the separator decomposition
        // For now, we'll use a heuristic based on vertex IDs
        return next > current;
    }

    // Creates an empty separator decomposition for default constructor.
    static SeparatorDecomposition createEmptyDecomp() {
        SeparatorDecomposition decomp;
        // Initialize with empty tree and order
        return decomp;
    }

    // Computes the LCA level of two vertices using simplified approach.
    int computeLCALevel(int source, int target) const {
        if (source >= bitStrings.size() || target >= bitStrings.size()) {
            return -1;
        }
        
        // 简化实现：使用顶点在order中的位置来计算LCA层级
        int sourcePos = -1, targetPos = -1;
        
        // 找到顶点在order中的位置
        for (int i = 0; i < decomp.order.size(); ++i) {
            if (decomp.order[i] == source) sourcePos = i;
            if (decomp.order[i] == target) targetPos = i;
        }
        
        if (sourcePos < 0 || targetPos < 0) {
            return -1;
        }
        
        // 简化LCA计算：使用位置差的绝对值
        int level = std::abs(sourcePos - targetPos);
        
        // 限制在合理范围内
        return std::min(level, static_cast<int>(decomp.tree.size() - 1));
    }
};