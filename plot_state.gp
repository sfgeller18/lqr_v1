# Set the title and labels
set title "Simulated State Trajectory"
set xlabel "Time Steps"
set ylabel "State Values"

# Specify the output format and file
set terminal png
set output "output/computedSimulatedStateTrajectory.png"

# Define the data file and column separation
data_file = "output/computedSimulatedStateTrajectory.csv"

# Define the plot style and legend entries
set style data lines

# Plot the data for each column
plot data_file using 1:2 with lines title "Column 1", \
     data_file using 1:3 with lines title "Column 2", \
     data_file using 1:4 with lines title "Column 3", \
     data_file using 1:5 with lines title "Column 4"
