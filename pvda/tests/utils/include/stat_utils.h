#ifndef STAT_UTILS_H
#define STAT_UTILS_H

// Chi-squared critical value for alpha = 0.05 depending on the degrees of freedom (values in source file)
extern const double chi_critical_05[100];

double jarque_bera(const double* x, int n);

#endif
