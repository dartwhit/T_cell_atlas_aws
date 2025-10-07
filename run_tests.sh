#!/bin/bash

# Test runner script for T-Cell Atlas AWS dynamic configuration
# Can be run locally or on EC2 instance

set -e  # Exit on any error

echo "🧪 T-Cell Atlas AWS - Dynamic Configuration Test Runner"
echo "======================================================"
echo "📅 Started at: $(date)"
echo "🖥️ Host: $(hostname)"
echo "📁 Working directory: $(pwd)"
echo ""

# Check if we're in the right directory
if [ ! -f "app_code/setup.R" ]; then
    echo "❌ Error: setup.R not found. Please run this script from the project root."
    exit 1
fi

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "❌ Error: R is not installed or not in PATH"
    exit 1
fi

echo "✅ R found: $(Rscript --version 2>&1 | head -1)"

# Install required packages if needed
echo ""
echo "📦 Checking R package dependencies..."
Rscript -e "
packages <- c('jsonlite', 'shiny', 'DT', 'ggplot2', 'dplyr', 'stringr', 'shinyjs', 'tidyr')
missing <- packages[!packages %in% installed.packages()[,'Package']]
if (length(missing) > 0) {
  cat('Installing missing packages:', paste(missing, collapse = ', '), '\n')
  install.packages(missing, repos = 'http://cran.us.r-project.org')
} else {
  cat('✅ All required packages are installed\n')
}
"

# Run the comprehensive test suite
echo ""
echo "🔬 Running comprehensive test suite..."
echo "--------------------------------------"

if Rscript tests/test_dataset_config.R; then
    echo ""
    echo "🎉 All tests passed!"
    
    # Show summary if test results file exists
    if [ -f "app_code/test_results.json" ]; then
        echo ""
        echo "📊 Test Summary:"
        if command -v jq &> /dev/null; then
            jq -r '.summary | "Total: \(.total), Passed: \(.passed), Failed: \(.failed), Success Rate: \(.success_rate * 100 | floor)%"' app_code/test_results.json
        else
            echo "   (Install jq for detailed summary)"
        fi
    fi
    
else
    echo ""
    echo "❌ Tests failed! Check output above for details."
    
    # Show failed tests if results file exists
    if [ -f "app_code/test_results.json" ] && command -v jq &> /dev/null; then
        echo ""
        echo "❌ Failed tests:"
        jq -r '.tests | to_entries[] | select(.value.passed == false) | "  • \(.key): \(.value.message)"' app_code/test_results.json
    fi
    
    exit 1
fi

# Run a basic smoke test of the application (if requested)
if [ "$1" = "--smoke-test" ]; then
    echo ""
    echo "💨 Running smoke test..."
    echo "----------------------"
    
    cd app_code
    
    # Test that the app can start without crashing
    timeout 30s Rscript -e "
        source('setup.R')
        source('server.R')
        source('ui.R')
        
        cat('✅ App components loaded successfully\n')
        cat('📊 Datasets available:', length(dataset_files), '\n')
        cat('🎯 Study levels configured for first dataset:\n')
        print(data_level_choices[[1]])
    " || echo "⚠️ Smoke test timed out or failed (this is expected without full data)"
    
    cd ..
fi

echo ""
echo "🏁 Test run completed at: $(date)"
echo "======================================================"