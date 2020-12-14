# Introduction 

This folder contains the source code (i.e., .py files) used in the paper “ALSO-X is Better Than CVaR: Convex Approximations for Chance Constrained
Programs Revisited.” The paper proposed two algorithms, termed "ALSO-X" and "ALSO-X+," to solve CCPs and DRCCPs, corresponding to subfolders Experiment 1-2 and subfolders Experiment 3-4. To successfully run the code in subfolders, you may need to install Python (i.e. .py files) and Matlab (i.e. .m files).
Below is the detailed explanation of each subfolder.


# Subfolder: Experiment 1
**Testing in CCPs.**
* **Purpose:**  Experiment 1 was implemented to test the performance of ALSO-X and ALSO-X+ for solving regular CCP. 
* **Source code info:** The code folder includes two .py files. "Linear CCP. py" and "Nonlinear CCP. py" are the codes for implementations in linear and nonlinear CCP, respectively.



# Subfolder: Experiment 2
**Comparisons in Covering CCPs.**
* **Purpose:** Experiment 2 was implemented to compare CVaR approximation, ALSO-X, ALSO-X+, and Relax-and-Scale Algorithm for solving covering CCPs.
* **Source code info:** The code folder includes one .py file. "Covering CCP. py" is the code for implementation.

# Subfolder: Experiment 3
**Testing in DRCCPs.**
* **Purpose:** Experiment 3 was implemented to test the performance of the worst-case ALSO-X and the worst-case ALSO-X+ for solving DRCCP. 
* **Source code info:** The code folder includes three .py files. "1 Wasserstein with 2 norm. py" and "1 Wasserstein with 2 norm & tuning. py" are the codes for implementations in DRCCP with 1-Wasserstein ambiguity set with 2-norm and the penalized version to improve the worst-case ALSO-X. "infinite Wasserstein with infinite norm. py" is the code for implementations in DRCCP with $`\infty`$-Wasserstein ambiguity set with $`\infty`$-norm. 


# Subfolder: Experiment 4
**Illustration Example in a DRCCP with Gaussian reference distribution.**
* **Purpose:** Experiment 4 was implemented to illustrate the comparisons in DRCCP with Gaussian reference distribution.
* **Source code info:** The code folder includes three .m files.  "eta_theta. m", "eta_epsilon. m", and "eta_comparison. m" are the codes for discussions in Figure4, Figure 5, and Figure 6, respectively.




Feel free to refer to the codes displayed in the subfolders. If used in your research, please cite our paper.

References: Jiang, N., Xie, W. (2020) ALSO-X is Better Than CVaR: Convex Approximations for Chance Constrained Programs Revisited. Available at Optimization Online.
http://www.optimization-online.org/DB_FILE/2020/12/8148.pdf
