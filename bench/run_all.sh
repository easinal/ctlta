#!/bin/bash
# Every number in note/CTNR_results_summary.md's main tables: 8 schemes x 6 graphs x 3 query
# loads, plus one customization run per graph and metric. Check the machine is idle
# (cat /proc/loadavg) first, or the numbers mean nothing.
#
#   ./run_all.sh [out dir]                defaults to out-all, then ./sum_all.py <dir>
#   SCHEMES="FHL CTNR" ./run_all.sh       only those schemes (CCH is always added: it is the
#                                         reference answer every other scheme is checked against)
#   GRAPHS="NY E" ./run_all.sh            only those graphs
#   NODE=0 THREADS=56 ./run_all.sh        one NUMA node only -- how every table in the doc so
#                                         far was measured (default: all threads, all nodes)
#   CACHE=0 ./run_all.sh                  do not cache the CCH
#   also overridable: EXE DD
#
# Output files are named   <graph>_<scheme>_<metric>_<load>.csv  (+ .log), e.g.
#   USA_FHL-inregion_time_rank.csv   = FHL + in-region hubs, travel-time metric, rank load
# with metric = time | dist   and   load = long (10000 uniform random pairs) | rank (Dijkstra
# rank buckets, 500 pairs each; the rank-256 bucket is the doc's "city trip" load).
# <graph>_customization_<metric>.csv holds the CCH and CTL customization times.
#
# Cost warning: the full matrix is 150 runs, a few hours. CTL on USA needs ~105 GB of RAM.
set -u
cd "$(dirname "$0")"
EXE=${EXE:-../Build/Release/Launchers/RunP2PAlgo}; DD=${DD:-/data/zwan018/dimacs}
OUT=${1:-out-all}; CACHE=${CACHE:-1}
mkdir -p "$OUT" query "$DD/cch"
[ -x $EXE ] || { echo "no $EXE -- build it: cmake --build ../Build/Release --target RunP2PAlgo"; exit 1; }

# Threads and memory placement. One node: threads and memory pinned together, so every memory
# access is local -- the setting the doc's tables use. All nodes: threads spread over the whole
# box, memory interleaved page by page across all nodes so no single node's bandwidth is the
# bottleneck; query latency then pays remote accesses for 3 of every 4 pages.
if [ -n "${NODE:-}" ]; then
  THREADS=${THREADS:-56}; NUMA=(numactl --cpunodebind=$NODE --membind=$NODE)
else
  THREADS=${THREADS:-$(nproc)}; NUMA=(numactl --interleave=all)
fi
# Record how this batch was run: the thread count and NUMA policy change the numbers, and a
# results table that does not say which setting produced it cannot be compared with anything.
{ echo "date:     $(date '+%Y-%m-%d %H:%M:%S')"
  echo "host:     $(hostname)"
  echo "threads:  $THREADS"
  echo "numa:     ${NUMA[*]}"
  echo "binary:   $(realpath $EXE)  md5 $(md5sum < $EXE | cut -c1-12)"
  echo "graphs:   ${GRAPHS:-NY FLA CAL E W USA}"
  echo "schemes:  ${SCHEMES:-all}"
  echo "loadavg at start: $(cut -d' ' -f1-3 /proc/loadavg)"
} > "$OUT/run-config.txt"
cat "$OUT/run-config.txt"

for g in ${GRAPHS:-NY FLA CAL E W USA}; do
  B=USA-road-d.$g

  # ---- Warm-up prefix: why we first run 2000 queries we then throw away ----------------
  # The rank query file is sorted by ascending Dijkstra rank, 500 queries per bucket. So the
  # "rank 256" bucket is exactly the first 500 queries the process ever runs -- caches, TLB
  # and page tables all cold. It carries the whole program's warm-up cost on its own: it
  # measures well above the steady state, and the SAME command re-run differs by 12% (USA)
  # to 39% (CAL).
  #
  # The fix is to run 2000 throwaway queries first to warm the index up. They are tagged
  # rank 0, and the summary script only reports the rank-256 bucket, so they drop out.
  #
  # The catch: the warm-up must use a DIFFERENT set of pairs (taken from the uniform-random
  # 10k file here). Replaying the rank queries themselves and keeping the second pass would
  # pull exactly the data those queries need into cache, and every bucket would come out
  # nearly twice as fast -- flattering, but not the real performance.
  RANKQ=query/$B.rank-with-warmup.csv
  if [ ! -s $RANKQ ]; then
    grep '^#' $DD/query/$B.rk.csv                       > $RANKQ   # comment header of the original
    echo origin,destination,dijkstra_rank              >> $RANKQ   # column header
    grep -vE '^#|^origin' $DD/query/$B.10k.csv | head -2000 |
      awk -F, '{print $1","$2",0"}'                    >> $RANKQ   # 2000 warm-up pairs, tagged rank 0
    grep -vE '^#|^origin' $DD/query/$B.rk.csv          >> $RANKQ   # the queries we actually measure
    echo "  wrote $RANKQ"
  fi

  # The CCH is metric-independent, so all runs of this graph share one cached copy (USA:
  # 18.5 s -> 2 s of start-up each, 3.1 GB). It lives next to the graph in $DD/cch. One file
  # per graph: the cache is not self-describing, so sharing a path across graphs is nonsense.
  CC=(); [ "$CACHE" = 1 ] && CC=(-cch-cache $DD/cch/$B.cch)

  # Per-graph knobs:
  #   CTNR             -ctnr-thresh auto: the cut is sized so the transit count T ~ 17.6 sqrt(n)
  #                    (the T x T table then costs the same per vertex on every graph; USA keeps
  #                    its 14 levels, FLA/CAL/E drop from 14 to 13)
  #   overlay levels   still 12 for NY, 14 for the rest
  #   subregion size   the in-region cut, roughly S/4 ~ S/6 of the `auto` region size
  TH=$([ $g = NY ] && echo 12 || echo 14)
  case $g in NY) SI=32;; FLA) SI=48;; CAL|E) SI=64;; W) SI=128;; USA) SI=256;; esac

  for load in "time_long $DD/query/$B.10k.csv" "time_rank $RANKQ" "dist_long $DD/query/$B.10k.csv -l"; do
    set -- $load; tag=$1 q=$2 metric=${3:-}   # third field is the optional -l (distance metric)
    for scheme in "CCH                    -a CCH-tree" \
                  "CTNR                   -a CTNR -ctnr-thresh auto -ctnr-prune-thresh 127" \
                  "FHL                    -a FHL -fhl-region-size auto" \
                  "FHL-inregion           -a FHL -fhl-region-size auto -fhl-subregion-size $SI" \
                  "FHL-toplabels          -a FHL -fhl-region-size auto -fhl-top-labels 3" \
                  "FHL-inregion-toplabels -a FHL -fhl-region-size auto -fhl-subregion-size $SI -fhl-top-labels 3" \
                  "FHL-overlay            -a FHL -ctnr-thresh $TH -fhl-overlay-every 3" \
                  "CTL                    -a CTL"; do
      set -- $scheme; name=$1; shift           # drop the name, pass the rest straight through
      case " CCH ${SCHEMES:-} " in *" $name "*) ;; *) [ -n "${SCHEMES:-}" ] && continue;; esac
      f=$OUT/${g}_${name}_$tag
      # The first run of each graph builds the CCH instead of loading it, so the report gets a
      # real construction time (loading the cache is ~3 s on USA, building it ~18.5 s).
      CCRUN=("${CC[@]}"); [ "$name" = CCH ] && [ "$tag" = time_long ] && CCRUN=()
      # Customization is reported from the time_long run: customize 3 times there and report the
      # median (-n, the same repetition format the *-custom modes use), for every scheme alike.
      REPS=(); [ "$tag" = time_long ] && REPS=(-n 3)
      CTNR_THREADS=$THREADS "${NUMA[@]}" $EXE "$@" $metric "${CCRUN[@]}" "${REPS[@]}" \
        -g $DD/bin/$B.gr.bin -s $DD/bisep/$B.strict_bisep.bin -d $q -o $f.csv > $f.log 2>&1
      echo "  $(basename $f)  EXIT=$?"         # no `set -e`: one bad run must not abort the batch
    done
  done

done
echo "done. summarise with: ./sum_all.py $OUT"
