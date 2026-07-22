Basin-of-attraction follow-up scripts
====================================

Put these files in the folder that contains your basin_of_attraction_outputs folder, then run in Julia:

    include("basin_clean_summary_figures.jl")
    include("basin_direction1_representative_trajectories.jl")
    include("basin_direction2_refined_tumor_nk_boundary.jl")

Required input for the first two scripts:
    basin_of_attraction_outputs/basin_all_runs.csv
    basin_of_attraction_outputs/basin_fraction_summary.csv

If Julia says the CSV cannot be found:
    1. Check your current directory in Julia:
          pwd()
    2. Change directory to the project folder:
          cd("/path/to/your/project/folder")
    3. Confirm the output folder exists:
          readdir()
          readdir("basin_of_attraction_outputs")

If Direction 2 is slow:
    Open basin_direction2_refined_tumor_nk_boundary.jl and change:
          GRID_N = 41
    to:
          GRID_N = 25
    for a faster test run.

Main outputs:
    basin_interpretation_figures/
    basin_direction1_trajectories/
    basin_direction2_refined_boundary/
