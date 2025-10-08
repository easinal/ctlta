#!/bin/bash

# CTL Workflow Script for DIMACS data
# Usage: ./run_p2p_workflow.sh [--build] [--rebuild] [--convert] [--bisep] [--test graph] [--check] [--force] [--graphs graph1,graph2,...]

# build from scratch
# sh run_p2p_workflow.sh --theta 5000 --threads 224 --cxx /usr/bin/g++-14 --cc /usr/bin/gcc --rebuild

# convert all graphs
# sh run_p2p_workflow.sh --convert

# convert specific graphs  
# sh run_p2p_workflow.sh --graphs USA-road-d.NY --convert

# manual convert example:
# ./Build/Release/RawData/ConvertGraph -s dimacs -d binary -c -scc -i /data/zwan018/dimacs/USA-road-d.NY -o /data/zwan018/dimacs/bin/USA-road-d.NY

# generate bisep files for all graphs
# sh run_p2p_workflow.sh --bisep

# generate bisep files for specific graphs
# sh run_p2p_workflow.sh --graphs USA-road-d.E,USA-road-d.W --bisep

# test single graph
# sh run_p2p_workflow.sh --test USA-road-d.E

# run all graphs
# sh run_p2p_workflow.sh

# run specific graphs
# sh run_p2p_workflow.sh --graphs USA-road-d.E,USA-road-d.W

# check results
# sh run_p2p_workflow.sh --check

# generate OD pairs only (no queries)
# sh run_p2p_workflow.sh --generate-od

# generate OD pairs with custom settings
# sh run_p2p_workflow.sh --od-pairs 5000 --od-seed 42 --od-lengths --generate-od

# force regenerate all files (even if they exist)
# sh run_p2p_workflow.sh --force
    
set -e

DIMACS_DIR="/data/zwan018/dimacs"
BIN_DIR="/data/zwan018/dimacs/bin"
BISEP_DIR="/data/zwan018/dimacs/bisep"
CTNR_DIR="/data/zwan018/dimacs/ctnr"
QUERY_DIR="/data/zwan018/dimacs/query"
RESULTS_DIR="/data/zwan018/dimacs/results"
EXE="Build/Release/Launchers/RunP2PAlgo"
CONVERT_EXE="Build/Release/RawData/ConvertGraph"
OD_EXE="Build/Release/RawData/GenerateODPairs"

# Default: run all graphs
SELECTED_GRAPHS=""
CTL_THETA=0
NUM_THREADS=224
CXX_COMPILER="/usr/bin/g++-14"
C_COMPILER="/usr/bin/gcc-14"
NUM_OD_PAIRS=100
OD_SEED=0
OD_USE_LENGTHS=false
FORCE_REGENERATE=false

log() { echo "[$(date '+%H:%M:%S')] $1"; }
success() { echo "✓ $1"; }
error() { echo "✗ $1"; }
warn() { echo "⚠ $1"; }

# Build project
build() {
    log "Building project with CTL_THETA=$CTL_THETA, NUM_THREADS=$NUM_THREADS..."
    [ ! -d "External/RoutingKit" ] && git submodule update --init
    [ ! -f "External/RoutingKit/lib/libroutingkit.so" ] && make -C External/RoutingKit lib/libroutingkit.so
    
    # Build cmake command
    local cmake_cmd="cmake -S . -B Build/Release -DCMAKE_BUILD_TYPE=Release -DCTL_THETA=$CTL_THETA -DNUM_THREADS=$NUM_THREADS"
    
    if [ -n "$CXX_COMPILER" ]; then
        cmake_cmd="$cmake_cmd -DCMAKE_CXX_COMPILER=$CXX_COMPILER"
        log "Using C++ compiler: $CXX_COMPILER"
    fi
    
    if [ -n "$C_COMPILER" ]; then
        cmake_cmd="$cmake_cmd -DCMAKE_C_COMPILER=$C_COMPILER"
        log "Using C compiler: $C_COMPILER"
    fi
    
    log "CMake command: $cmake_cmd"
    eval "$cmake_cmd"
    cmake --build Build/Release --target RunP2PAlgo ConvertGraph GenerateODPairs
    success "Build complete"
}

# Clean and rebuild project from scratch
rebuild() {
    log "Cleaning and rebuilding project from scratch..."
    
    # Clean build directory
    if [ -d "Build" ]; then
        log "Removing existing build directory..."
        rm -rf Build
    fi
    
    # Clean RoutingKit build files
    if [ -d "External/RoutingKit" ]; then
        log "Cleaning RoutingKit build files..."
        rm -rf External/RoutingKit/build External/RoutingKit/lib 2>/dev/null || true
    fi
    
    # Reinitialize submodules
    log "Reinitializing submodules..."
    git submodule update --init --recursive
    
    # Build RoutingKit from scratch
    log "Building RoutingKit from scratch..."
    make -C External/RoutingKit lib/libroutingkit.so
    
    # Configure and build project
    log "Configuring project with CTL_THETA=$CTL_THETA, NUM_THREADS=$NUM_THREADS..."
    
    # Build cmake command
    local cmake_cmd="cmake -S . -B Build/Release -DCMAKE_BUILD_TYPE=Release -DCTL_THETA=$CTL_THETA -DNUM_THREADS=$NUM_THREADS"
    
    if [ -n "$CXX_COMPILER" ]; then
        cmake_cmd="$cmake_cmd -DCMAKE_CXX_COMPILER=$CXX_COMPILER"
        log "Using C++ compiler: $CXX_COMPILER"
    fi
    
    if [ -n "$C_COMPILER" ]; then
        cmake_cmd="$cmake_cmd -DCMAKE_C_COMPILER=$C_COMPILER"
        log "Using C compiler: $C_COMPILER"
    fi
    
    log "CMake command: $cmake_cmd"
    eval "$cmake_cmd"
    
    log "Building executables..."
    cmake --build Build/Release --target RunP2PAlgo ConvertGraph GenerateODPairs
    
    success "Rebuild complete!"
}

# Convert DIMACS graphs to binary format
convert_graphs() {
    log "Converting DIMACS graphs to binary format..."
    mkdir -p "$BIN_DIR"
    
    if [ -n "$SELECTED_GRAPHS" ]; then
        # Convert selected graphs only
        for graph_name in $(echo "$SELECTED_GRAPHS" | tr ',' ' '); do
            local bin_file="$BIN_DIR/${graph_name}.gr.bin"
            
            if [ ! -f "$bin_file" ] || [ "$FORCE_REGENERATE" = true ]; then
                if [ "$FORCE_REGENERATE" = true ] && [ -f "$bin_file" ]; then
                    log "Force regenerating $graph_name..."
                    rm -f "$bin_file"
                else
                    log "Converting $graph_name..."
                fi
                $CONVERT_EXE -s dimacs -d binary -c -scc -a length travel_time lat_lng -i "$DIMACS_DIR/$graph_name" -o "$BIN_DIR/$graph_name"
            else
                log "Binary file already exists: $bin_file"
            fi
        done
    else
        # Convert all graphs
        for dimacs_file in "$DIMACS_DIR"/*.dist.gr; do
            [ -f "$dimacs_file" ] || continue
            local base=$(basename "$dimacs_file" .dist.gr)
            local bin_file="$BIN_DIR/${base}.gr.bin"
            
            if [ ! -f "$bin_file" ] || [ "$FORCE_REGENERATE" = true ]; then
                if [ "$FORCE_REGENERATE" = true ] && [ -f "$bin_file" ]; then
                    log "Force regenerating $base..."
                    rm -f "$bin_file"
                else
                    log "Converting $base..."
                fi
                $CONVERT_EXE -s dimacs -d binary -c -scc -a length travel_time lat_lng -i "$DIMACS_DIR/$base" -o "$BIN_DIR/$base"
            else
                log "Binary file already exists: $bin_file"
            fi
        done
    fi
    success "Graph conversion complete"
}

# Generate bisep files for all graphs
generate_all_bisep() {
    log "Generating bisep files..."
    mkdir -p "$BISEP_DIR"
    
    if [ -n "$SELECTED_GRAPHS" ]; then
        # Generate bisep for selected graphs only
        for graph_name in $(echo "$SELECTED_GRAPHS" | tr ',' ' '); do
            local graph="$BIN_DIR/${graph_name}.gr.bin"
            local bisep_file="$BISEP_DIR/${graph_name}.strict_bisep.bin"
            
            if [ -f "$graph" ]; then
                if [ ! -f "$bisep_file" ] || [ "$FORCE_REGENERATE" = true ]; then
                    if [ "$FORCE_REGENERATE" = true ] && [ -f "$bisep_file" ]; then
                        log "Force regenerating bisep for $graph_name..."
                        rm -f "$bisep_file"
                    else
                        log "Generating bisep for $graph_name..."
                    fi
                    $EXE -a CTL -g "$graph" -o "$bisep_file" -b 30
                else
                    log "Bisep file already exists: $bisep_file"
                fi
            else
                warn "Binary graph file not found: $graph"
            fi
        done
    else
        # Generate bisep for all graphs
        for graph in "$BIN_DIR"/*; do
            [ -f "$graph" ] || continue
            local base=$(basename "$graph" .gr.bin)
            local bisep_file="$BISEP_DIR/${base}.strict_bisep.bin"
            
            if [ ! -f "$bisep_file" ] || [ "$FORCE_REGENERATE" = true ]; then
                if [ "$FORCE_REGENERATE" = true ] && [ -f "$bisep_file" ]; then
                    log "Force regenerating bisep for $base..."
                    rm -f "$bisep_file"
                else
                    log "Generating bisep for $base..."
                fi
                $EXE -a CTL -g "$graph" -o "$bisep_file" -b 30
            else
                log "Bisep file already exists: $bisep_file"
            fi
        done
    fi
    success "Bisep generation complete"
}

# Generate CTNR files for all graphs
generate_all_ctnr() {
    log "Generating CTNR files..."
    mkdir -p "$CTNR_DIR"
    
    if [ -n "$SELECTED_GRAPHS" ]; then
        # Generate CTNR for selected graphs only
        for graph_name in $(echo "$SELECTED_GRAPHS" | tr ',' ' '); do
            local graph="$BIN_DIR/${graph_name}.gr.bin"
            local ctnr_file="$CTNR_DIR/${graph_name}.ctnr.bin"
            
            if [ -f "$graph" ]; then
                if [ ! -f "$ctnr_file" ] || [ "$FORCE_REGENERATE" = true ]; then
                    if [ "$FORCE_REGENERATE" = true ] && [ -f "$ctnr_file" ]; then
                        log "Force regenerating CTNR for $graph_name..."
                        rm -f "$ctnr_file"
                    else
                        log "Generating CTNR for $graph_name..."
                    fi
                    $EXE -a CTNR -g "$graph" -o "$ctnr_file"
                else
                    log "CTNR file already exists: $ctnr_file"
                fi
            else
                warn "Binary graph file not found: $graph"
            fi
        done
    else
        # Generate CTNR for all graphs
        for graph in "$BIN_DIR"/*; do
            [ -f "$graph" ] || continue
            local base=$(basename "$graph" .gr.bin)
            local ctnr_file="$CTNR_DIR/${base}.ctnr.bin"
            
            if [ ! -f "$ctnr_file" ] || [ "$FORCE_REGENERATE" = true ]; then
                if [ "$FORCE_REGENERATE" = true ] && [ -f "$ctnr_file" ]; then
                    log "Force regenerating CTNR for $base..."
                    rm -f "$ctnr_file"
                else
                    log "Generating CTNR for $base..."
                fi
                $EXE -a CTNR -g "$graph" -o "$ctnr_file"
            else
                log "CTNR file already exists: $ctnr_file"
            fi
        done
    fi
    success "CTNR generation complete"
}

# Generate bisep files if not exist
generate_bisep() {
    local graph="$1"
    local base=$(basename "$graph" .gr.bin)
    local bisep_file="$BISEP_DIR/${base}.strict_bisep.bin"
    
    if [ ! -f "$bisep_file" ] || [ "$FORCE_REGENERATE" = true ]; then
        if [ "$FORCE_REGENERATE" = true ] && [ -f "$bisep_file" ]; then
            log "Force regenerating bisep file for $base..."
            rm -f "$bisep_file"
        else
            log "Generating bisep file for $base..."
        fi
        mkdir -p "$BISEP_DIR"
        $EXE -a CTL -g "$graph" -o "$bisep_file" -b 30
    else
        log "Bisep file already exists: $bisep_file"
    fi
}

# Generate OD pairs if not exist
generate_od_pairs() {
    local graph="$1"
    local base=$(basename "$graph" .gr.bin)
    local query_file="$QUERY_DIR/${base}.csv"
    
    if [ ! -f "$query_file" ] || [ "$FORCE_REGENERATE" = true ]; then
        if [ "$FORCE_REGENERATE" = true ] && [ -f "$query_file" ]; then
            log "Force regenerating $NUM_OD_PAIRS OD pairs for $base..."
            rm -f "$query_file"
        else
            log "Generating $NUM_OD_PAIRS OD pairs for $base..."
        fi
        mkdir -p "$QUERY_DIR"
        
        local cmd="$OD_EXE -g $graph -o $query_file -n $NUM_OD_PAIRS -s $OD_SEED"
        if [ "$OD_USE_LENGTHS" = true ]; then
            cmd="$cmd -l"
        fi
        
        log "OD generation command: $cmd"
        eval "$cmd"
    else
        log "Query file already exists: $query_file"
    fi
}

# Generate OD pairs only (no queries)
generate_od_pairs_only() {
    log "Generating OD pairs for all graphs..."
    
    if [ -n "$SELECTED_GRAPHS" ]; then
        # Generate OD pairs for selected graphs only
        for graph_name in $(echo "$SELECTED_GRAPHS" | tr ',' ' '); do
            local graph="$BIN_DIR/${graph_name}.gr.bin"
            if [ -f "$graph" ]; then
                generate_od_pairs "$graph"
            else
                warn "Binary graph file not found: $graph"
            fi
        done
    else
        # Generate OD pairs for all graphs
        for graph in "$BIN_DIR"/*; do
            [ -f "$graph" ] || continue
            generate_od_pairs "$graph"
        done
    fi
    
    success "OD pairs generation complete!"
}

# Process single graph
process_graph() {
    local graph="$1"
    local base=$(basename "$graph" .gr.bin)
    local bisep_file="$BISEP_DIR/${base}.strict_bisep.bin"
    local ctnr_file="$CTNR_DIR/${base}.ctnr.bin"
    local query_file="$QUERY_DIR/${base}.csv"
    local result_file="$RESULTS_DIR/${base}_ctl_results"
    
    log "Processing $base with CTL..."
    
    # Generate bisep file if needed
    generate_bisep "$graph"
    
    # Generate CTNR file if needed
    if [ ! -f "$ctnr_file" ] || [ "$FORCE_REGENERATE" = true ]; then
        if [ "$FORCE_REGENERATE" = true ] && [ -f "$ctnr_file" ]; then
            log "Force regenerating CTNR for $base..."
            rm -f "$ctnr_file"
        else
            log "Generating CTNR for $base..."
        fi
        mkdir -p "$CTNR_DIR"
        $EXE -a CTNR -g "$graph" -o "$ctnr_file"
    else
        log "CTNR file already exists: $ctnr_file"
    fi
    
    # Generate OD pairs if needed
    generate_od_pairs "$graph"
    
    # Check file sizes and validity
    log "Checking files..."
    log "Graph: $graph ($(du -h "$graph" | cut -f1))"
    log "Bisep: $bisep_file ($(du -h "$bisep_file" | cut -f1))"
    log "CTNR: $ctnr_file ($(du -h "$ctnr_file" | cut -f1))"
    log "Query: $query_file ($(wc -l < "$query_file") lines)"
    
    # Create results directory
    mkdir -p "$RESULTS_DIR"
    
    # Run CTL queries with error handling
    log "Running CTL queries..."
    if $EXE -a CTL -g "$graph" -s "$bisep_file" -d "$query_file" -o "$result_file"; then
        success "CTL queries completed for $base"
    else
        error "CTL queries failed for $base (exit code: $?)"
        return 1
    fi
    
    # Run CTL customization with error handling (commented out for performance)
    # log "Running CTL customization..."
    # if $EXE -a CTL-custom -g "$graph" -s "$bisep_file" -o "${result_file}_custom.csv" -n 10; then
    #     success "CTL customization completed for $base"
    # else
    #     error "CTL customization failed for $base (exit code: $?)"
    #     return 1
    # fi
    
    success "Completed $base"
}

# Run Dijkstra baseline for comparison
run_dijkstra_baseline() {
    local graph="$1"
    local base=$(basename "$graph" .gr.bin)
    local query_file="$QUERY_DIR/${base}.csv"
    local dijkstra_file="$RESULTS_DIR/${base}_dijkstra_baseline.csv"
    
    if [ ! -f "$query_file" ]; then
        warn "Query file not found: $query_file"
        return 1
    fi
    
    log "Running Dijkstra baseline for $base..."
    if $EXE -a Dij -g "$graph" -d "$query_file" -o "$dijkstra_file"; then
        success "Dijkstra baseline completed for $base"
        return 0
    else
        error "Dijkstra baseline failed for $base"
        return 1
    fi
}

# Compare CTL results with Dijkstra baseline
compare_results() {
    local ctl_file="$1"
    local dijkstra_file="$2"
    local base="$3"
    
    if [ ! -f "$ctl_file" ] || [ ! -f "$dijkstra_file" ]; then
        echo "$base | N/A | N/A | Missing files"
        return 1
    fi
    
    # Compare distances (first column)
    local mismatches=$(paste <(tail -n +2 "$ctl_file" | cut -d',' -f1) <(tail -n +2 "$dijkstra_file" | cut -d',' -f1) | awk -F'\t' '$1 != $2' | wc -l)
    local total_queries=$(tail -n +2 "$ctl_file" | wc -l)
    
    if [ "$mismatches" -eq 0 ]; then
        echo "$base | $total_queries | 0 | Perfect match"
        return 0
    else
        local error_rate=$((mismatches * 100 / total_queries))
        echo "$base | $total_queries | $mismatches ($error_rate%) | Distance mismatch"
        return 1
    fi
}

# Check results and generate statistics
check_results() {
    log "Checking results and generating statistics..."
    
    local total_queries=0
    local total_time=0
    local failed_queries=0
    local processed_graphs=0
    local accuracy_checked=0
    local perfect_matches=0
    
    echo "=========================================="
    echo "CTL Results Summary"
    echo "=========================================="
    echo "Graph Name | Queries | Avg Time (ns) | Status"
    echo "------------------------------------------"
    
    for result_file in "$RESULTS_DIR"/*_ctl_results.csv; do
        [ -f "$result_file" ] || continue
        
        local base=$(basename "$result_file" _ctl_results.csv)
        local query_file="$QUERY_DIR/${base}.csv"
        local graph_file="$BIN_DIR/${base}.gr.bin"
        
        if [ ! -f "$query_file" ]; then
            echo "$base | N/A | N/A | Missing query file"
            continue
        fi
        
        # Count queries
        local num_queries=$(tail -n +2 "$result_file" | wc -l)
        local num_expected=$(tail -n +2 "$query_file" | wc -l)
        
        if [ "$num_queries" -ne "$num_expected" ]; then
            echo "$base | $num_queries/$num_expected | N/A | Query count mismatch"
            ((failed_queries++))
            continue
        fi
        
        # Calculate average time
        local avg_time=0
        if [ "$num_queries" -gt 0 ]; then
            avg_time=$(tail -n +2 "$result_file" | cut -d',' -f2 | awk '{sum+=$1; count++} END {if(count>0) print sum/count; else print 0}')
        fi
        
        # Check for errors (negative distances or times)
        local errors=$(tail -n +2 "$result_file" | awk -F',' '$1<0 || $2<0' | wc -l)
        
        if [ "$errors" -gt 0 ]; then
            echo "$base | $num_queries | $avg_time | $errors errors found"
            ((failed_queries++))
        else
            echo "$base | $num_queries | $avg_time | OK"
            ((total_queries += num_queries))
            ((total_time += avg_time * num_queries))
        fi
        
        ((processed_graphs++))
    done
    
    echo "------------------------------------------"
    echo "Total processed graphs: $processed_graphs"
    echo "Total queries: $total_queries"
    echo "Failed graphs: $failed_queries"
    
    if [ "$total_queries" -gt 0 ]; then
        local overall_avg=$((total_time / total_queries))
        echo "Overall average query time: $overall_avg ns"
    fi
    
    # Accuracy check with Dijkstra baseline
    echo ""
    echo "=========================================="
    echo "Accuracy Check (CTL vs Dijkstra)"
    echo "=========================================="
    echo "Graph Name | Queries | Mismatches | Status"
    echo "------------------------------------------"
    
    for result_file in "$RESULTS_DIR"/*_ctl_results.csv; do
        [ -f "$result_file" ] || continue
        
        local base=$(basename "$result_file" _ctl_results.csv)
        local dijkstra_file="$RESULTS_DIR/${base}_dijkstra_baseline.csv"
        local graph_file="$BIN_DIR/${base}.gr.bin"
        
        # Run Dijkstra baseline if not exists
        if [ ! -f "$dijkstra_file" ] && [ -f "$graph_file" ]; then
            log "Running Dijkstra baseline for $base..."
            run_dijkstra_baseline "$graph_file"
        fi
        
        # Compare results
        if [ -f "$dijkstra_file" ]; then
            if compare_results "$result_file" "$dijkstra_file" "$base"; then
                ((perfect_matches++))
            fi
            ((accuracy_checked++))
        else
            echo "$base | N/A | N/A | No Dijkstra baseline"
        fi
    done
    
    echo "------------------------------------------"
    echo "Accuracy checked graphs: $accuracy_checked"
    echo "Perfect matches: $perfect_matches"
    
    if [ "$accuracy_checked" -gt 0 ]; then
        local accuracy_rate=$((perfect_matches * 100 / accuracy_checked))
        echo "Accuracy rate: $accuracy_rate%"
    fi
    
    # Check customization results
    echo ""
    echo "=========================================="
    echo "Customization Results"
    echo "=========================================="
    echo "Graph Name | Custom Runs | Avg Time (μs)"
    echo "------------------------------------------"
    
    for custom_file in "$RESULTS_DIR"/*_ctl_results_custom.csv; do
        [ -f "$custom_file" ] || continue
        
        local base=$(basename "$custom_file" _ctl_results_custom.csv)
        local num_runs=$(tail -n +2 "$custom_file" | wc -l)
        
        if [ "$num_runs" -gt 0 ]; then
            local avg_custom=$(tail -n +2 "$custom_file" | cut -d',' -f3 | awk '{sum+=$1; count++} END {if(count>0) print sum/count; else print 0}')
            echo "$base | $num_runs | $avg_custom"
        else
            echo "$base | 0 | N/A"
        fi
    done
    
    echo "=========================================="
    
    if [ "$failed_queries" -eq 0 ] && [ "$perfect_matches" -eq "$accuracy_checked" ]; then
        success "All results are correct and match Dijkstra baseline!"
    elif [ "$failed_queries" -eq 0 ]; then
        warn "Results generated but some don't match Dijkstra baseline"
    else
        warn "Found $failed_queries graphs with issues"
    fi
}

# Get list of graphs to process
get_graphs() {
    if [ -n "$SELECTED_GRAPHS" ]; then
        # Use selected graphs
        for graph_name in $(echo "$SELECTED_GRAPHS" | tr ',' ' '); do
            local graph_file="$BIN_DIR/${graph_name}.gr.bin"
            if [ -f "$graph_file" ]; then
                echo "$graph_file"
            else
                warn "Graph file not found: $graph_file"
            fi
        done
    else
        # Use all graphs
        for graph in "$BIN_DIR"/*; do
            [ -f "$graph" ] && echo "$graph"
        done
    fi
}

# Main execution
main() {
    if [ "$1" = "--build" ]; then
        build
        exit 0
    fi
    
    if [ "$1" = "--convert" ]; then
        [ ! -f "$CONVERT_EXE" ] && { error "Build first: $0 --build"; exit 1; }
        convert_graphs
        exit 0
    fi
    
    [ ! -f "$EXE" ] && { error "Build first: $0 --build"; exit 1; }
    
    local graphs_to_process=($(get_graphs))
    if [ ${#graphs_to_process[@]} -eq 0 ]; then
        error "No graphs found to process"
        exit 1
    fi
    
    log "Processing ${#graphs_to_process[@]} graphs..."
    for graph in "${graphs_to_process[@]}"; do
        process_graph "$graph"
    done
    
    success "All graphs processed! Results in $RESULTS_DIR"
}

# Parse command line arguments
ACTION=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --build)
            ACTION="build"
            shift
            ;;
        --rebuild)
            ACTION="rebuild"
            shift
            ;;
        --convert)
            ACTION="convert"
            shift
            ;;
        --bisep)
            ACTION="bisep"
            shift
            ;;
        --ctnr)
            ACTION="ctnr"
            shift
            ;;
        --test)
            ACTION="test"
            TEST_GRAPH="$2"
            shift 2
            ;;
        --check)
            ACTION="check"
            shift
            ;;
        --generate-od)
            ACTION="generate-od"
            shift
            ;;
        --theta)
            CTL_THETA="$2"
            log "CTL_THETA set to: $CTL_THETA"
            shift 2
            ;;
        --threads)
            NUM_THREADS="$2"
            log "NUM_THREADS set to: $NUM_THREADS"
            shift 2
            ;;
        --cxx)
            CXX_COMPILER="$2"
            log "C++ compiler set to: $CXX_COMPILER"
            shift 2
            ;;
        --cc)
            C_COMPILER="$2"
            log "C compiler set to: $C_COMPILER"
            shift 2
            ;;
        --od-pairs)
            NUM_OD_PAIRS="$2"
            log "Number of OD pairs set to: $NUM_OD_PAIRS"
            shift 2
            ;;
        --od-seed)
            OD_SEED="$2"
            log "OD generation seed set to: $OD_SEED"
            shift 2
            ;;
        --od-lengths)
            OD_USE_LENGTHS=true
            log "OD generation will use physical lengths"
            shift
            ;;
        --force)
            FORCE_REGENERATE=true
            log "Force regenerate all files (even if they exist)"
            shift
            ;;
        --graphs)
            SELECTED_GRAPHS="$2"
            shift 2
            ;;
        *)
            error "Unknown option: $1"
            echo "Usage: $0 [--build] [--rebuild] [--convert] [--bisep] [--ctnr] [--test graph] [--check] [--generate-od] [--force] [--theta value] [--threads num] [--cxx compiler] [--cc compiler] [--od-pairs num] [--od-seed seed] [--od-lengths] [--graphs graph1,graph2,...]"
            echo "Examples:"
            echo "  $0 --build"
            echo "  $0 --rebuild  # Clean rebuild from scratch"
            echo "  $0 --theta 10000 --threads 224 --cxx /usr/bin/g++-14 --cc /usr/bin/gcc --rebuild"
            echo "  $0 --convert"
            echo "  $0 --bisep"
            echo "  $0 --ctnr"
            echo "  $0 --generate-od  # Generate OD pairs only (no queries)"
            echo "  $0 --od-pairs 5000 --od-seed 42 --generate-od  # Generate 5000 OD pairs with seed 42"
            echo "  $0 --force  # Force regenerate all files (even if they exist)"
            echo "  $0 --test USA-road-d.E  # Test single graph"
            echo "  $0 --check  # Check results and show statistics"
            echo "  $0  # Run all graphs"
            echo "  $0 --graphs USA-road-d.E,USA-road-d.W  # Run specific graphs (auto-complete to .gr.bin)"
            exit 1
            ;;
    esac
done

# Execute action after all parameters are parsed
case $ACTION in
    "build")
        build
        exit 0
        ;;
    "rebuild")
        rebuild
        exit 0
        ;;
    "convert")
        [ ! -f "$CONVERT_EXE" ] && { error "Build first: $0 --build"; exit 1; }
        convert_graphs
        exit 0
        ;;
    "bisep")
        [ ! -f "$EXE" ] && { error "Build first: $0 --build"; exit 1; }
        generate_all_bisep
        exit 0
        ;;
    "ctnr")
        [ ! -f "$EXE" ] && { error "Build first: $0 --build"; exit 1; }
        generate_all_ctnr
        exit 0
        ;;
    "test")
        [ ! -f "$EXE" ] && { error "Build first: $0 --build"; exit 1; }
        test_graph="$BIN_DIR/$TEST_GRAPH.gr.bin"
        if [ ! -f "$test_graph" ]; then
            error "Test graph not found: $test_graph"
            exit 1
        fi
        log "Testing single graph: $TEST_GRAPH"
        process_graph "$test_graph"
        exit 0
        ;;
    "check")
        check_results
        exit 0
        ;;
    "generate-od")
        [ ! -f "$OD_EXE" ] && { error "Build first: $0 --build"; exit 1; }
        generate_od_pairs_only
        exit 0
        ;;
    "")
        # No action specified, run main workflow
        main
        ;;
    *)
        error "Unknown action: $ACTION"
        exit 1
        ;;
esac
