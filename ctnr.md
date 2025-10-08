# Ideas for Customizable Transit Node Routing

## 1. Choosing Transit Nodes
* original paper: top k ranked vertices in CH
* here: top k levels of SD (level = distance from root of SD)

## 2. Finding Access Nodes of a Single Vertex v (probably not a great idea, see next section for better one)
* explore upward search space of v until you hit transit node and prune
* all transit nodes hit are access nodes
* for metric-dependent:
	* after CCH perfect customization -> will allow smaller sets of access nodes due to removed shortcuts
	* same thing as metric-independent, except Dijkstra-based instead of BFS

## 3. Finding Access Nodes of All Vertices Systematically
* let A(v) be set of access nodes for vertex v
* initialize A(v) = {v} for every transit node
* then, top-down sweep starting at level k+1
* assume we know access nodes of all vertices up to and including level l
* for vertex v in level l+1:
	* A(v) := union of A(w) for every upward neighbor
* for metric-dependent (i.e. during customization and possibly after CCH perfect customization):
	* all transit nodes have distance 0 to themselves
	* assume two upward neighbors x,y of v both have entry for a transit node u with distances d(x,u), d(y,u)
	* then, d(v,u) := min { d(v,x) + d(x,u), d(v,y) + d(y,u) }
	* (so union but take minimum of distances if more than one neighbor has an entry)

## 4. Locality Filter
* store location of each vertex in tree decomposition as bit string -> e.g., 0 for left child, 1 for right child at every level
* then, an XOR operation is enough to find level of lowest common ancestor (LCA) in tree decomposition
* if level(LCA(s,t)) > k (i.e., further from the root than any transit node), then query is local, otherwise non-local
	* for non-local queries, every path will contain at least one transit node, so it's okay to only use transit nodes for them
	* for local queries, the shortest path may contain a transit node but does not have to so both TN query and local query are needed

## 5. Local Query
* local queries are likely to have shortest path that does not contain transit node (but it may)
* run graph search (e.g., regular CCH query) AND run TN query 
* local graph search can be pruned to only consider vertices that are less important than transit nodes since that part will be covered by TN query
	* Dijkstra-based: prune whenever a transit node is hit
	* elimination-tree-based: stop when first transit node is hit (all remaining vertices will be at least as high in SD, so will also be transit nodes)
* could think about a multi-level idea where the local graph search is replaced by another level of TNR or something -> multi-level TNR

## 6. Peter's Idea for Optimization of Non-Local TN Queries
* let l: = level(LCA(s,t))
* if l <= k, then the path from s to t will contain a transit node
* however, the path will only contain transit nodes whose level is at most l (since the separator at l separates s and t)
* so TN query for s,t only has to consider any access node a of s and t for which level(a) <= l (instead of a <= k)
* let A_l(s) = { a in A(s) | level(a) <= l } and A_l(t) = analogous
* without the optimization, we need |A(s)| \* |A(t)| reads from the distance table
* with the optimization it's |A_l(s)| \* |A_l(t)| reads which might be much smaller


## 7. Overview of Three Phases
* Metric-independent Preprocessing:
	* build separator decomposition (e.g. using InertialFlowCutter), store bit strings for location of each vertex
	* decide level k based on number of targeted transit nodes
	* compute CCH shortcut edges (metric-independent)
* Customization:
	* run CCH basic customization (computes weights of shortcuts)
	* run CCH perfect customization (removes shortcuts that are not needed for the given metric)
	* compute distance table of transit nodes (in parallel)
	* compute access nodes along with distances systematically using approach from point 3. above
* Query:
	* check if query is local using locality filter (very fast)
	* if query is local: 
		* run pruned local graph search
		* then run full TN query 
			* (can we do something better here? Peter's idea does not work for this case. Maybe prune with result from graph search, e.g. ignore access nodes with larger distances or something?)
	* if query is non-local:
		* run TN query with Peter's pruning


## Comparison to CTL
+ space usage (depending on theta and number of transit nodes)
+ query time for unimportant vertices
	* query time not decided by importance of vertex but by local/non-local
- query time in general
	* three-hop instead of two-hop needs number of reads from distance table quadratic in number of access nodes
		* could be reduced through LCA-based pruning of TN queries
	* perhaps there are many access nodes compared to regular TNR because of metric-independent shortcuts
		* could be okay with perfect customization (i.e., access nodes decided upon during customization after CCH perfect customization)
- customization time in general?
	* might actually be okay compared to labeling with correct systematic approach
	* should generally be much less work for TNR rather than HL