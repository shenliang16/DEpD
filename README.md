# DEpD
Code of paper "Robust Line Segment Mismatch Removal Using Point-pair Representation and Gaussian-Uniform Mixture Formulation"

Mistake correction in the paper:

Eq. 26

{\rm{M_2}} = \frac{1}{{25}}\sqrt {\sum\nolimits_{i = 1}^{25} {\left\| {\rho \left( {{\bf{\hat Hy}}_i^*} \right) - \rho \left( {{\bf{Hy}}_i^*} \right)} \right\|_2^2} } 

is not correct, which should be

{\rm{M_2}} = \sqrt {\frac{1}{{25}}\sum\nolimits_{i = 1}^{25} {\left\| {\rho \left( {{\bf{\hat Hy}}_i^*} \right) - \rho \left( {{\bf{Hy}}_i^*} \right)} \right\|_2^2} }

We are sorry for the mistake.
