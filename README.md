# uncertain-control-co-design-SASA-example
This code creates the reults for the deterministic control co-design, stochastic in expectation uncertain control co-design, and worst-case robust uncertain control co-design for the simple SASA case study.
The results are associated with the following article: 

-Azad,S., and Herber, D. R., 2022."Investigations into uncertain control co-design implementations for stochastic in expectation and worst-case robust". In International Mechanical Engineering Congress & Exposition, no. IMECE2022-95229

Start by openning the Main.m files which contains nine different case studies associated with the porblems discussed in the above-mentioned article.

- case 0 solves the deterministic simple SASA problem

- case 1 solves the Stochastic in expectation simple SASA problem in nested architecture using a multiple control approach through Monte Carlo Simulation

- case 2 solves the Stochastic in expectation simple SASA problem in nested architecture using a multiple control approach through generalized Polynomial Chaos

- case 3 solves the worst case robust simple SASA problem using direct single shooting and a single control approach through constraint relaxation

- case 4 solves the worst case robust simple SASA problem using direct single shooting and a single control approach with penalty terms

- case 5 solves the worst case robust simple SASA problem using a multiple control approach through through polytopic uncertainties

- case 6 and 7 solve case 5 for various sizes of uncertainties 

- case 7 solves closed-loop system response

- case 8 implements a robust multi-stage model predicitve control approach using polytopic uncertainties
  
# External Dependencies
* [MATLAB Figure Workflow](https://github.com/danielrherber/matlab-figure-workflow)
* [DTQP Project](https://github.com/danielrherber/dt-qp-project)
* [Export_Fig](https://github.com/altmany/export_fig)
* [mfoldername](https://www.mathworks.com/matlabcentral/fileexchange/40397-mfoldername)
* [Quadrature - forked from Aero-matlab](https://github.com/wme7/Aero-matlab)
