ctnr_access_node_pruning.md - Markdown
# Customizing Access Nodes with Pruning
* basic approach:
	* customize distance table
	* initialize access node sets for transit nodes: every node is its own single access node with distance 0
	* propagate access nodes top-down

* notation:
	* d(x,y) means overall shortest-path distance from x to y
	* dup(x,y) means upward distance from x to y
	* ddown(x,y) means downward distance from x to y

## Propagation Process
* assume we know access nodes of all vertices up to and including level l
* for vertex v in level l+1:
	* A(v) := union of A(w) for every upward neighbor, resolve duplicate occurrences by taking minimum distance
		* assume two upward neighbors x,y of v both have access node w with distances dup(x,w), dup(y,w)
		* then, dup(v,w) := min { dup(v,x) + dup(x,w), dup(v,y) + dup(y,w) }
* new step: prune
	* let A(v) = {w_1,...,w_r}, ordered by increasing rank
	* right after constructing A(v) during the sweep, do the following
	* for i = 2,...,r do:
		* for j = 1,...,i-1 do:
			* if dup(v,w_j) + d(w_j,w_i) <= dup(v,w_i), then w_i is dominated by w_j and we can remove w_i from A(v)
				* d(w_j,w_i) can be read from distance table

## Correctness:
* assume w_j dominates w_i for A(v) and a shortest path from v to u uses access node w_i for v and access node x for u, so:
	* d(v,u) = dup(v,w_i) + d(w_i,x) + ddown(x,u)

* transform this based on the domination:		
	* d(v,u) = dup(v,w_i) + d(w_i,x) 			  + ddown(x,u)
	        >= dup(v,w_j) + d(w_j,w_i) + d(w_i,x) + ddown(x,u)		// use domination
			>= dup(v,w_j) + d(w_j,x) 			  + ddown(x,u)		// d(w_j,x) <= d(w_j,w_i) + d(w_i,x)
			>= d(v,u) 												// any path from v to u is at least as long as shortest path

* thus, a path using access node w_j for v and access node x for u is also a shortest path
* this path will be found by the query since w_j remains an access node of v
	* (because w_j is checked for pruning before w_i so if w_j dominates w_i, w_j is sure to stay in A(v))