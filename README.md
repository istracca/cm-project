READ ME

Folder:
- graph_instances: This folder contains various DIMACS files (graph instances) used to generate different matrices for our problem. 

Functions:
- create_matrix_from_dimacs.m: This function reads a graph instance from a DIMACS format file and creates the corresponding matrix representations.
- lanczos.m: This function performs the Lanczos process to generate an orthogonal basis and tridiagonal matrix for a given matrix.
- qr_iterative.m: This function performs iterative QR decomposition to update the QR factors of a matrix incrementally.
- trim_schur_component.m: This function trims the Schur component matrix by either performing strict trimming (removing all non-diagonal elements) or soft trimming (keeping the k largest diagonal elements above a threshold).
- custom_minres: This function implements a custom MINRES algorithm.
- custom_minres_preconditioned.m: This function implements a custom MINRES algorithm with preconditioning to solve linear systems.

Scripts:
- custom_vs_matlab.m: This script compares the performance of custom_minres and MATLAB's MINRES by measuring computation time, residuals, and iterations for various graph instances.
- custom_vs_matlab_plot.m: This script compares the performance of custom_minres and MATLAB's MINRES by plotting the evolution of relative residuals for both methods.
- preconditioned_vs_original_plot.m: This script compares the performance of preconditioned and original custom_minres by plotting the evolution of relative residuals for both methods.
- preconditioned_vs_original.m: This script compares the performance of preconditioned and original custom_minres by measuring computation time, residuals, and iterations for various graph instances.
