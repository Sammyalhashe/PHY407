# PHY407, Fall 2015, Lab 4, Q1b
# Author: DUONG, BANG CHI

# Test the accuracy and timing of 3 methods: GaussElim, PartialPivot, and LU decomposition : Solve Ax=v

# Import modules
from numpy import linspace, empty, dot, mean
from numpy.random import rand
from numpy.linalg import solve
from time import clock
from SolveLinear import GaussElim, PartialPivot
from matplotlib.pyplot import figure, plot, title, xlabel, ylabel, legend, yscale, show

# Create an N_array
N_array = linspace(5,500)

# Define time arrays for 3 methods
time_GaussElim = empty(len(N_array))
time_PartialPivot = empty(len(N_array))
time_LU = empty(len(N_array))
time_index = 0

# Define error arrays for 3 methods
error_GaussElim = empty(len(N_array))
error_PartialPivot = empty(len(N_array))
error_LU = empty(len(N_array))
error_index = 0

for N_index in N_array:
    # Set A and v randomly from N_array
    v = rand(N_index)
    A = rand(N_index, N_index)

    # -------- For Gaussian elimination-----------------------
    start = clock()
    x_GaussElim = GaussElim(A,v)
    end = clock()
    # Add elements to the time array
    time_GaussElim[time_index] = end - start
    # Check if the answer is correct
    v_sol_GaussElim = dot(A, x_GaussElim)
    # Add elements to the error array
    error_GaussElim[error_index] = mean(abs(v-v_sol_GaussElim))

     # -------- For Partial Pivot--------------------------------
    start = clock()
    x_PartialPivot = PartialPivot(A,v)
    end = clock()
    # Add elements to the time array
    time_PartialPivot[time_index] = end - start
    # Check if the answer is correct
    v_sol_PartialPivot = dot(A, x_PartialPivot)
    # Add elements to the error array
    error_PartialPivot[error_index] = mean(abs(v-v_sol_PartialPivot))

    # -------- For LU decomposition--------------------------
    start = clock()
    x_LU = solve(A,v)
    end = clock()
    # Add elements to the time array
    time_LU[time_index] = end - start
    # Check if the answer is correct
    v_sol_LU = dot(A, x_LU)
    # Add elements to the error array
    error_LU[error_index] = mean(abs(v-v_sol_LU))

    time_index += 1
    error_index += 1

# Plot for time taken to operate each method
figure()
plot(N_array, time_GaussElim, label = 'Gaussian Elimination')
plot(N_array, time_PartialPivot, label = 'Partial Pivoting')
plot(N_array, time_LU, label = 'LU Decompostion')
legend(loc = "lower right")
title('Running time for three methods')
xlabel('N')
ylabel('Time taken in log scale')
yscale('log')


# Plot for error in each method
figure()
plot(N_array, error_GaussElim, label = 'Gaussian Elimination')
plot(N_array, error_PartialPivot, label = 'Partial Pivoting')
plot(N_array, error_LU, label = 'LU Decompostion')
legend(loc = "lower right")
title('Error for three methods')
xlabel('N')
ylabel('Error in log scale')
yscale('log')

show()

