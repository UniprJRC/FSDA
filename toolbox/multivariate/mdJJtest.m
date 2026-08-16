function out = mdJJtest(Y, varargin)
%mdJJtest Jamshidian-Jalal test for MCAR in incomplete continuous data.
%
%<a href="matlab: docsearchFS('mdJJtest')">Link to the help function</a>
%
%  mdJJtest implements the diagnostic framework proposed by Jamshidian and
%  Jalal (2010). The Hawkins, Neyman and k-sample Anderson-Darling
%  calculations follow the modern reimplementation in function mcar of the
%  R package mice, which is itself based on TestMCARNormality from the
%  MissMech package. The present function is not a line-by-line reproduction
%  of MissMech: it uses FSDA joint-normal stochastic imputation, supports
%  multiple user-supplied completed data sets, reports median results across
%  imputations and includes additional numerical regularization.
%
%
%  Required input arguments:
%
%    Y :           Input data. Matrix or table. n x p data matrix possibly
%                  containing missing values (NaN's). Rows of Y represent
%                  observations and columns represent variables.
%                  Data Types - single | double | table
%
%  Optional input arguments:
%
%    imputed :     Completed data sets. Empty value, cell array or 3D array.
%                  If supplied, imputed can be:
%                  (i) an n x p x M numeric array; or
%                  (ii) a cell array containing M n x p numeric matrices or
%                  tables. Observed entries should not be changed. If empty,
%                  mdEM and mdImputeStochastic are used internally.
%                  The default value is [].
%                  Example - 'imputed',Yimp
%                  Data Types - double | single | cell
%
% nimputations :   Number of stochastic imputations. Positive integer.
%                  This option is used only when imputed is empty.
%                  The default value is 5.
%                  Example - 'nimputations',10
%                  Data Types - double | single
%
%       minn :     Minimum retained pattern size. Integer greater than or
%                  equal to 2. Missingness patterns containing fewer than
%                  minn observations are removed. Therefore, every retained
%                  pattern contains at least minn observations.
%                  The default value is 6.
%                  Example - 'minn',10
%                  Data Types - double | single
%
%     method :     Test method. Character vector or string.
%                  Possible values are:
%                  'auto'          : Hawkins' test is computed first. If at
%                                    least one Hawkins p-value is smaller
%                                    than alpha, the nonparametric
%                                    Anderson-Darling test is also computed;
%                  'hawkins'       : only Hawkins' test is reported;
%                  'nonparametric' : the Anderson-Darling rank test is
%                                    always reported.
%                  The default value is 'auto'.
%                  Example - 'method','nonparametric'
%                  Data Types - char | string
%
%     nsimul :     Number of Monte Carlo replications used to calibrate the
%                  fourth-order Neyman smooth test in small patterns.
%                  The default value is 10000.
%                  Example - 'nsimul',50000
%                  Data Types - double | single
%
%  usechisq :      Threshold for the asymptotic Neyman reference law.
%                  For a pattern containing at least usechisq observations,
%                  the chi-square approximation is used instead of Monte
%                  Carlo calibration. The default value is 30.
%                  Example - 'usechisq',50
%                  Data Types - double | single
%
%      alpha :     Significance level. Scalar in the interval (0,1).
%                  The default value is 0.05.
%                  Example - 'alpha',0.01
%                  Data Types - double | single
%
%
%    maxiter :     Maximum number of iterations for mdEM. Positive integer.
%                  This option is used only when imputed is empty.
%                  The default value is 200.
%                  Example - 'maxiter',500
%                  Data Types - double | single
%
%        tol :     Convergence tolerance for mdEM. Positive scalar.
%                  This option is used only when imputed is empty.
%                  The default value is 1e-7.
%                  Example - 'tol',1e-8
%                  Data Types - double | single
%
%      ridge :     Numerical regularization. Nonnegative scalar.
%                  A small diagonal ridge is used to stabilize the pooled
%                  covariance matrix in Hawkins' test. The default value is
%                  1e-10.
%                  Example - 'ridge',1e-8
%                  Data Types - double | single
%
%        msg :     Display test results. Boolean.
%                  If true, a concise summary is displayed.
%                  The default value is true.
%                  Example - 'msg',false
%                  Data Types - logical | double
%
%      plots :     Produce graphical summaries of the diagnostic results.
%                  Boolean. The default value is false.
%
%                  If plots is true, the function creates:
%
%                  1) A figure containing bar charts of the available
%                     imputation-specific p-values.
%
%                     If Hawkins' test is computed, the first bar chart
%                     displays the elements of
%                     $\mathtt{out.hawkinsP}$. Bar $m$ is the Hawkins
%                     p-value obtained from completed data set $m$. A
%                     horizontal reference line is drawn at alpha. The
%                     title reports the median Hawkins p-value, which is
%                     the value used as the main Hawkins result, together
%                     with the percentage of imputation-specific p-values
%                     smaller than alpha.
%
%                     If the Anderson-Darling test is computed, a second
%                     bar chart displays the elements of
%                     $\mathtt{out.ADpvalue}$. Bar $m$ is the
%                     Anderson-Darling p-value obtained from completed data
%                     set $m$. A horizontal reference line is drawn at
%                     alpha. The title reports the median Anderson-Darling
%                     p-value, which is used as the main result when the
%                     Anderson-Darling branch is selected, together with
%                     the percentage of imputation-specific p-values
%                     smaller than alpha.
%
%                     When both tests are available, the two bar charts are
%                     shown in separate panels of the same figure. The
%                     horizontal axis identifies the completed data sets,
%                     while the vertical axis shows the corresponding
%                     p-values.
%
%                     Under method='auto', the Anderson-Darling test is
%                     activated when at least one Hawkins p-value is smaller
%                     than alpha. This activation criterion is distinct
%                     from the final decision, which is based on the median
%                     Anderson-Darling p-value when that branch is computed.
%
%                  2) A graphical representation of the missingness
%                     patterns retained for the MCAR test. Rows correspond
%                     to retained missingness patterns and columns to
%                     variables.
%
%                     The graphical representation is selected
%                     automatically by mdpattern. For a moderate number of
%                     patterns, a balloon plot is used: large circles denote
%                     missing entries and small filled dots denote observed
%                     entries. When the number of patterns is large, a
%                     heatmap is used, in which red cells denote missing
%                     entries and blue cells denote observed entries.
%
%                     The left axis reports the frequency of each displayed
%                     pattern, while the right axis reports the number of
%                     missing variables in that pattern. Variable names are
%                     shown on the top axis and the number of missing
%                     entries for each variable is shown on the bottom
%                     axis.
%
%                     Observations belonging to patterns containing fewer
%                     than minn observations are excluded before the graph
%                     is constructed. Consequently, the displayed
%                     frequencies and missing-value totals refer only to
%                     observations retained for the MCAR test.
%
%                  3) A Pareto plot of the mean contribution of each
%                     retained missingness pattern to every available test
%                     statistic. Contributions are first computed
%                     separately for the M completed data sets and then
%                     averaged pattern by pattern.
%
%                     Bars are sorted from the largest to the smallest mean
%                     contribution. The line on the right axis reports the
%                     cumulative percentage of the total mean statistic,
%                     and a dashed horizontal line marks 80 percent. Labels
%                     such as P3 (n=12) identify the original pattern number
%                     and its frequency. These identifiers are the same as
%                     those used in the missingness-pattern figure.
%
%                     Hawkins and Anderson-Darling contributions are shown
%                     in separate panels when both tests are available.
%                     The plot is diagnostic: it reveals whether the global
%                     statistic is driven mainly by a few patterns or is
%                     distributed across many patterns.
%
%                  Setting plots to false suppresses all graphical output
%                  but does not affect the numerical results returned in
%                  out.
%
%                  Example - 'plots',true
%                  Data Types - logical | double
%
%  Output:
%
%    out :         Structure containing the following fields:
%
%          out.stat          = main test statistic selected by the chosen
%                              method.
%
%                              If method='hawkins', or if method='auto' and
%                              the Anderson-Darling branch is not activated,
%                              out.stat is the median of the
%                              imputation-specific combined Hawkins
%                              statistics stored in
%                              $\mathtt{out.hawkinsStat}$:
%
%                              \[
%                                  \mathtt{out.stat}
%                                  =
%                                  \mathrm{median}
%                                  \left\{
%                                  H^{(1)},\ldots,H^{(M)}
%                                  \right\}.
%                              \]
%
%                              Here, $H^{(m)}$ is Fisher's combined Hawkins
%                              statistic for completed data set $m$.
%
%                              If method='nonparametric', or if method='auto'
%                              and the Anderson-Darling branch is activated,
%                              out.stat is the median of the
%                              imputation-specific Anderson-Darling
%                              statistics stored in $\mathtt{out.ADstat}$:
%
%                              \[
%                                  \mathtt{out.stat}
%                                  =
%                                  \mathrm{median}
%                                  \left\{
%                                  A^{2(1)}_{g,n},\ldots,
%                                  A^{2(M)}_{g,n}
%                                  \right\}.
%                              \]
%
%                              Thus, out.stat is a descriptive median across
%                              completed data sets; it is not obtained using
%                              Rubin's rules or another formal
%                              multiple-imputation pooling procedure.
%
%          out.pvalue        = main p-value associated with out.stat.
%
%                              If the Hawkins branch supplies the main
%                              result, out.pvalue is the median of the
%                              imputation-specific Hawkins p-values stored
%                              in $\mathtt{out.hawkinsP}$:
%
%                              \[
%                                  \mathtt{out.pvalue}
%                                  =
%                                  \mathrm{median}
%                                  \left\{
%                                  p_{\mathrm{H}}^{(1)},\ldots,
%                                  p_{\mathrm{H}}^{(M)}
%                                  \right\}.
%                              \]
%
%                              If the Anderson-Darling branch supplies the
%                              main result, out.pvalue is the median of the
%                              imputation-specific Anderson-Darling p-values
%                              stored in $\mathtt{out.ADpvalue}$:
%
%                              \[
%                                  \mathtt{out.pvalue}
%                                  =
%                                  \mathrm{median}
%                                  \left\{
%                                  p_{\mathrm{AD}}^{(1)},\ldots,
%                                  p_{\mathrm{AD}}^{(M)}
%                                  \right\}.
%                              \]
%
%                              Under method='auto', the Anderson-Darling
%                              branch is activated when at least one element
%                              of $\mathtt{out.hawkinsP}$ is smaller than
%                              alpha. Once this branch is activated,
%                              out.stat and out.pvalue are based on the
%                              Anderson-Darling results, not on the Hawkins
%                              medians.
%
%                              The comparison of out.pvalue with alpha gives
%                              the reported decision. The complete vectors
%                              $\mathtt{out.hawkinsP}$ and
%                              $\mathtt{out.ADpvalue}$ should also be
%                              inspected because they show the variability
%                              of the results across completed data sets.
%
%                              The reported median p-value is a descriptive
%                              summary and is not a formally pooled
%                              multiple-imputation p-value.
%          out.df            = degrees of freedom of the main test. This is
%                              empty for the Anderson-Darling test;
%          out.method        = method used;
%          out.alpha         = significance level;
%          out.hawkinsStat   = M x 1 vector containing the combined Hawkins
%                              statistics, one for each imputation;
%          out.hawkinsDF     = degrees of freedom of the combined Hawkins
%                              statistic;
%          out.hawkinsP      = M x 1 vector of Hawkins p-values;
%          out.ADstat        = M x 1 vector of Anderson-Darling statistics,
%                              when computed;
%          out.ADpvalue      = M x 1 vector of Anderson-Darling p-values,
%                              when computed;
%          out.meanHawkinsContribution = g x 1 vector containing the mean
%                              contribution of each retained missingness
%                              pattern to the Hawkins statistic across the M
%                              completed data sets. Element r is the mean of
%                              $-2\log(p_r^{(m)})$ over m. The field is empty
%                              when Hawkins' test is not computed;
%          out.meanADContribution = g x 1 vector containing the mean
%                              contribution of each retained missingness
%                              pattern to the Anderson-Darling statistic
%                              across the M completed data sets. The field is
%                              empty when the Anderson-Darling test is not
%                              computed;
%          out.patterns      = g x p logical matrix containing the retained
%                              missingness patterns. A true entry denotes a
%                              missing variable;
%          out.patternCounts = g x 1 vector containing the frequencies of
%                              the retained missingness patterns;
%          out.patternInfo   = table containing pattern frequency and number
%                              of missing variables. Row names identify the
%                              missingness patterns;
%          out.npatterns     = number of retained missingness patterns;
%          out.idxPatterns   = nused x 1 pattern membership vector;
%          out.removedRows   = n x 1 logical vector identifying observations
%                              removed because of sparse patterns;
%          out.removedPatterns = matrix containing removed patterns;
%          out.loc           = EM estimate of location used for the default
%                              stochastic imputations. Empty if imputed is
%                              supplied;
%          out.cov           = EM estimate of covariance used for the default
%                              stochastic imputations. Empty if imputed is
%                              supplied;
%          out.details       = cell array with detailed results for each
%                              completed data set. When Hawkins' test is
%                              computed,
%                              $\mathtt{out.details\{m\}.HawkinsContributions}$
%                              contains the g pattern contributions for
%                              completed data set m. When the
%                              Anderson-Darling test is computed, the
%                              corresponding values are stored in
%                              $\mathtt{out.details\{m\}.AndersonDarling.}$
%                              $\mathtt{GroupContributions}$;
%          out.interpretation = concise interpretation of the test;
%          out.minn          = minimum retained pattern size. Missingness
%                              patterns containing fewer than
%                              $\mathtt{out.minn}$ observations are removed
%                              before the test is computed;
%
%          out.nsimul        = number of Monte Carlo replications used to
%                              approximate the null distribution of the
%                              fourth-order Neyman smooth statistic for
%                              patterns containing fewer than
%                              $\mathtt{out.usechisq}$ observations;
%
%          out.usechisq      = pattern-size threshold for the Neyman
%                              uniformity test. For a pattern containing at
%                              least $\mathtt{out.usechisq}$ observations,
%                              the asymptotic chi-square distribution with
%                              four degrees of freedom is used. For smaller
%                              patterns, Monte Carlo calibration based on
%                              $\mathtt{out.nsimul}$ replications is used;
%
%          out.nimputations  = number $M$ of completed data sets analysed.
%                              This is equal to the number of stochastic
%                              imputations generated internally when
%                              imputed is empty, or to the number of
%                              completed data sets supplied through
%                              imputed;
%
%          out.n             = number of observations in the original input
%                              data set, before removing observations
%                              belonging to sparse missingness patterns;
%
%          out.nused         = number of observations retained and used in
%                              the test after removing observations
%                              belonging to missingness patterns containing
%                              fewer than $\mathtt{out.minn}$
%                              observations;
%
%          out.p             = number of variables in the input data set;
%
%
%  More About:
%
%  The Jamshidian-Jalal procedure assesses the MCAR assumption by comparing
%  observations belonging to different missingness patterns. Because the
%  test requires complete vectors, it is applied separately to each
%  completed data set. The missingness-pattern groups are always defined
%  using the original incomplete data.
%
%  Let a retained completed data set contain $n$ observations, $p$
%  variables and $g$ missingness patterns. For pattern $r$, let $n_r$ be its
%  frequency, let $Y_{ri}$ denote completed observation $i$ in that pattern,
%  and let
%
%  \[
%      \overline{Y}_r
%      =
%      \frac{1}{n_r}\sum_{i=1}^{n_r}Y_{ri}
%  \]
%
%  be the corresponding pattern-specific mean.
%
%  Hawkins' transformation (Hawkins, 1981) begins with the pooled
%  within-pattern covariance matrix
%
%  \[
%      \widehat{\Sigma}_{P}
%      =
%      \frac{1}{n-g}
%      \sum_{r=1}^{g}
%      \sum_{i=1}^{n_r}
%      (Y_{ri}-\overline{Y}_r)
%      (Y_{ri}-\overline{Y}_r)^{\mathsf T}.
%  \]
%
%  Thus, variation between the pattern-specific means is removed before the
%  common covariance matrix is estimated. Under MCAR and multivariate
%  normality, the covariance matrices of the pattern groups should be
%  homogeneous.
%
%  For observation $i$ in pattern $r$, the squared Mahalanobis distance from
%  the pattern mean is
%
%  \[
%      d_{ri}
%      =
%      (Y_{ri}-\overline{Y}_r)^{\mathsf T}
%      \widehat{\Sigma}_{P}^{-1}
%      (Y_{ri}-\overline{Y}_r).
%  \]
%
%  Define
%
%  \[
%      h_{ri}=n_r d_{ri}.
%  \]
%
%  Hawkins' transformation is then
%
%  \[
%      F_{ri}
%      =
%      \frac{(n-g-p)h_{ri}}
%      {p\left\{(n_r-1)(n-g)-h_{ri}\right\}}.
%  \]
%
%  Under multivariate normality and equality of the covariance matrices
%  across missingness patterns, $F_{ri}$ has an $F$ distribution with $p$
%  and $n-g-p$ degrees of freedom. Consequently, the upper-tail probability
%
%  \[
%      A_{ri}
%      =
%      \Pr\left\{
%      F_{p,n-g-p}\geq F_{ri}
%      \right\}
%  \]
%
%  should follow a uniform distribution on $(0,1)$ within every missingness
%  pattern.
%
%
%  Uniformity is assessed separately in each pattern using a fourth-order
%  Neyman smooth test (Neyman, 1937; Rayner and Best, 1989). Let
%  $\phi_1,\ldots,\phi_4$ denote the first four orthonormal shifted
%  Legendre polynomials. The fourth order means that the first four
%  components of the orthonormal polynomial expansion are retained. These
%  components detect smooth departures from uniformity in four different
%  directions. For pattern $r$, the statistic is
%
%  \[
%      N_r
%      =
%      \frac{1}{n_r}
%      \sum_{\ell=1}^{4}
%      \left\{
%      \sum_{i=1}^{n_r}
%      \phi_{\ell}(A_{ri})
%      \right\}^{2}.
%  \]
%
%  For a pattern containing at least  $\mathtt{usechisq}$ observations,
%  $N_r$ is compared with a chi-square distribution having four degrees of
%  freedom:
%
%  \[
%      N_r \mathrel{\dot{\sim}} \chi^2_4.
%  \]
%
%  For smaller patterns, the reference distribution is approximated by
%  generating nsimul samples of size $n_r$ from the uniform distribution
%  and recomputing the Neyman statistic for each sample.
%
%  Let $p_r$ be the Neyman smooth-test $p$-value obtained for pattern $r$.
%  The pattern-specific results are combined using Fisher's method
%  (Fisher, 1932), giving the statistic
%
%
%  \[
%      H
%      =
%      -2\sum_{r=1}^{g}\log(p_r).
%  \]
%
%  Under the joint null hypothesis,
%
%  \[
%      H \mathrel{\dot{\sim}} \chi^2_{2g}.
%  \]
%
%  A small combined Hawkins $p$-value indicates that at least one pattern
%  departs from the expected uniform behaviour. Such a rejection can be
%  caused either by covariance heterogeneity, which provides evidence
%  against MCAR under the assumptions of the procedure, or by departure
%  from multivariate normality.
%
%  When method='hawkins', only this combined Hawkins test is reported. When
%  method='auto', Hawkins' test is computed first. If at least one
%  imputation-specific Hawkins p-value is smaller than alpha, the function
%  also applies the k-sample Anderson-Darling test. When
%  method='nonparametric', the Anderson-Darling test is always computed.
%
%  The k-sample Anderson-Darling test of Scholz and Stephens (1987)
%  compares the distributions of the transformed values $F_{ri}$ across
%  the $g$ retained missingness patterns. It extends the weighted
%  empirical-distribution criterion introduced by Anderson and Darling
%  (1952) to the k-sample problem. Its null hypothesis is
%
%  \[
%      H_0:\mathcal{F}_1=\cdots=\mathcal{F}_g,
%  \]
%
%  where $\mathcal{F}_r$ denotes the distribution of $F_{ri}$ in pattern
%  $r$. Under MCAR, these pattern-specific distributions should be equal.
%
%  Let
%
%  \[
%      z_1<\cdots<z_L
%  \]
%
%  denote the distinct ordered values in the pooled sample of transformed
%  observations, excluding the largest pooled value. For each $z_l$, let
%  $h_l$ be its multiplicity in the pooled sample, let $H_l$ be the number
%  of pooled observations not greater than $z_l$, and let $M_{rl}$ be the
%  number of observations from pattern $r$ not greater than $z_l$.
%
%  The contribution of pattern $r$ to the Anderson-Darling statistic is
%
%  \[
%      A_r
%      =
%      \frac{1}{n_r}
%      \sum_{l=1}^{L}
%      h_l
%      \frac{
%      \left(nM_{rl}-n_rH_l\right)^2
%      }{
%      H_l(n-H_l)
%      }.
%  \]
%
%  The complete k-sample statistic is
%
%  \[
%      A^2_{g,n}
%      =
%      \frac{1}{n}
%      \sum_{r=1}^{g}A_r.
%  \]
%
%  The denominator $H_l(n-H_l)$ gives relatively high weight to
%  distributional differences occurring in the tails. The multiplicities
%  $h_l$ provide an adjustment for tied transformed values.
%
%  Under the null hypothesis, the expected value of $A^2_{g,n}$ is
%  approximately $g-1$. The function computes its finite-sample variance,
%  denoted by $\sigma^2_{g,n}$, and forms the standardized statistic
%
%  \[
%      T_{g,n}
%      =
%      \frac{
%      A^2_{g,n}-(g-1)
%      }{
%      \sqrt{\sigma^2_{g,n}}
%      }.
%  \]
%
%  Large positive values of $T_{g,n}$ indicate stronger differences among
%  the pattern-specific distributions and therefore correspond to small
%  Anderson-Darling $p$-values.
%
%  Following the approximation of Scholz and Stephens (1987), the
%  Anderson-Darling $p$-value is not computed from a chi-square or an $F$
%  reference distribution. The implementation uses critical values
%  corresponding to the five upper-tail probabilities
%
%  \[
%      p_l
%      \in
%      \{0.25,0.10,0.05,0.025,0.01\}.
%  \]
%
%  For each probability $p_l$, the corresponding critical value is
%  approximated as a function of the number of retained patterns:
%
%  \[
%      q_l(g)
%      =
%      b_{0l}
%      +
%      \frac{b_{1l}}{\sqrt{g-1}}
%      +
%      \frac{b_{2l}}{g-1}.
%  \]
%
%  The tabulated probabilities are transformed to log-odds:
%
%  \[
%      c_l
%      =
%      \log\left(
%      \frac{1-p_l}{p_l}
%      \right).
%  \]
%
%  A cubic spline interpolates the value $c$ corresponding to the observed
%  standardized statistic $T_{g,n}$. The approximate upper-tail probability
%  is then
%
%  \[
%      p_{\mathrm{AD}}
%      =
%      \frac{1}{1+\exp(c)}.
%  \]
%
%  Therefore, $\mathtt{out.ADpvalue}(m)$ is the approximate
%  Anderson-Darling p-value obtained from completed data set $m$. A small
%  value indicates that the distributions of the Hawkins-transformed
%  quantities differ across the retained missingness patterns.
%
%  Because the calculation uses spline interpolation, and may use spline
%  extrapolation when $T_{g,n}$ lies outside the tabulated range, very small
%  or very large Anderson-Darling p-values should be interpreted as
%  approximate tail probabilities.
%
%  The k-sample Anderson-Darling test is distribution free and is used to
%  help distinguish a Hawkins rejection caused mainly by nonnormality from a
%  rejection caused by differences among the missingness-pattern
%  distributions. In the automatic procedure, a significant
%  Anderson-Darling result provides evidence against MCAR within the
%  MCAR-versus-MAR framework. A nonsignificant Anderson-Darling result
%  following a significant Hawkins test suggests that nonnormality is the
%  more plausible explanation.
%
%  When $M$ completed data sets are analysed, the entire procedure is
%  repeated separately for each imputation. The imputation-specific
%  statistics and p-values are retained in  $\mathtt{out.hawkinsStat}$,
%   $\mathtt{out.hawkinsP}$,  $\mathtt{out.ADstat}$ and
%   $\mathtt{out.ADpvalue}$.
%
%  The main reported statistic and $p$-value are the medians of the
%  imputation-specific quantities associated with the selected branch:
%
%  \[
%      T_{\mathrm{reported}}
%      =
%      \mathop{\mathrm{median}}_{m=1,\ldots, M}
%      T^{(m)},
%      \qquad
%      p_{\mathrm{reported}}
%      =
%      \mathop{\mathrm{median}}_{m=1,\ldots, M}
%      p^{(m)}.
%  \]
%
%  In particular, when the Anderson-Darling branch supplies the main
%  result,
%
%  \[
%      \mathtt{out.pvalue}
%      =
%      \mathop{\mathrm{median}}
%      \left\{
%      p_{\mathrm{AD}}^{(1)},\ldots,
%      p_{\mathrm{AD}}^{(M)}
%      \right\}.
%  \]
%
%  This median is a descriptive summary across imputations. It is not a
%  Rubin-rules or Meng-Rubin pooled p-value. The complete vectors of
%  imputation-specific p-values should therefore also be examined.
%
%  The global Hawkins and Anderson-Darling statistics can also be
%  decomposed into nonnegative pattern-specific contributions. For
%  Hawkins' test, let $p_r^{(m)}$ be the Neyman smooth-test p-value for
%  retained pattern r in completed data set m. Its contribution is
%
%  \[
%      C_{\mathrm{H},r}^{(m)}
%      =
%      -2\log\left(p_r^{(m)}\right),
%  \]
%
%  and therefore
%
%  \[
%      H^{(m)}
%      =
%      \sum_{r=1}^{g}C_{\mathrm{H},r}^{(m)}.
%  \]
%
%  The mean Hawkins contribution returned for pattern r is
%
%  \[
%      \overline{C}_{\mathrm{H},r}
%      =
%      \frac{1}{M}
%      \sum_{m=1}^{M}
%      \left[-2\log\left(p_r^{(m)}\right)\right].
%  \]
%
%  For the Anderson-Darling test, let
%  $C_{\mathrm{AD},r}^{(m)}$ denote the normalized group contribution
%  returned for pattern r in completed data set m. These contributions
%  satisfy
%
%  \[
%      A_{g,n}^{2(m)}
%      =
%      \sum_{r=1}^{g}C_{\mathrm{AD},r}^{(m)},
%  \]
%
%  and their mean is
%
%  \[
%      \overline{C}_{\mathrm{AD},r}
%      =
%      \frac{1}{M}
%      \sum_{m=1}^{M}C_{\mathrm{AD},r}^{(m)}.
%  \]
%
%  Consequently, the sums of the returned mean contributions satisfy
%
%  \[
%      \sum_{r=1}^{g}\overline{C}_{\mathrm{H},r}
%      =
%      \frac{1}{M}\sum_{m=1}^{M}H^{(m)},
%      \qquad
%      \sum_{r=1}^{g}\overline{C}_{\mathrm{AD},r}
%      =
%      \frac{1}{M}\sum_{m=1}^{M}A_{g,n}^{2(m)}.
%  \]
%
%  Thus, the mean contributions decompose the mean imputation-specific
%  statistic. They do not generally decompose $\mathtt{out.stat}$, because
%  the latter is the median of the imputation-specific statistics.
%
%  To construct the Pareto plot, let
%
%  \[
%      \overline{C}_{(1)}
%      \geq
%      \overline{C}_{(2)}
%      \geq
%      \cdots
%      \geq
%      \overline{C}_{(g)}
%  \]
%
%  denote the mean contributions sorted from largest to smallest. The
%  cumulative percentage after the first r sorted patterns is
%
%  \[
%      100
%      \frac{\sum_{j=1}^{r}\overline{C}_{(j)}}
%      {\sum_{j=1}^{g}\overline{C}_{(j)}}.
%  \]
%
%  Large bars identify patterns that contribute most strongly, on average
%  across completed data sets, to the corresponding global statistic. A
%  large contribution is a diagnostic indication and should not be
%  interpreted as a separate formal rejection for that pattern.
%
%  Option nsimul affects only the Monte Carlo calibration of the
%  pattern-specific Neyman tests in the Hawkins branch. It does not affect
%  the Anderson-Darling p-values.
%
%  When imputed is empty, mdJJtest estimates location and covariance using
%  mdEM and generates stochastic completed data sets using
%  mdImputeStochastic. This is a joint multivariate-normal imputation and is
%  not identical to the default chained normal-regression imputation used
%  by mice::mcar. To compare the MATLAB and R implementations using exactly
%  the same completed data sets, supply them through option imputed.
%
%  The procedure addresses the MCAR-versus-MAR framework. It cannot
%  distinguish MCAR from MNAR. Consequently, a nonsignificant result does
%  not establish MCAR and does not rule out an MNAR mechanism.
%
%  See also: mdEM, mdImputeStochastic, mdLittletest, mdMCARtest
%
%  References:
%
%  Anderson, T. W. and Darling, D. A. (1952), "Asymptotic Theory of
%  Certain Goodness-of-Fit Criteria Based on Stochastic Processes",
%  The Annals of Mathematical Statistics, Vol. 23, No. 2, pp. 193-212.
%  DOI: 10.1214/aoms/1177729437.
%
%  Fisher, R. A. (1932), "Statistical Methods for Research Workers",
%  4th ed., Oliver and Boyd, Edinburgh.
%
%  Hawkins, D. M. (1981), "A New Test for Multivariate Normality and
%  Homoscedasticity", Technometrics, Vol. 23, No. 1, pp. 105-110.
%  DOI: 10.1080/00401706.1981.10486244.
%
%  Jamshidian, M. and Jalal, S. (2010), "Tests of Homoscedasticity,
%  Normality, and Missing Completely at Random for Incomplete Multivariate
%  Data", Psychometrika, Vol. 75, pp. 649-674.
%
%  Jamshidian, M., Jalal, S. and Jansen, C. (2014), "MissMech: An R Package
%  for Testing Homoscedasticity, Multivariate Normality, and Missing
%  Completely at Random (MCAR)", Journal of Statistical Software, Vol. 56,
%  No. 6, pp. 1-31. DOI: 10.18637/jss.v056.i06.
%
%  Neyman, J. (1937), "'Smooth' Test for Goodness of Fit",
%  Skandinavisk Aktuarietidskrift, Vol. 20, pp. 149-199.
%
%  Rayner, J. C. W. and Best, D. J. (1989), "Smooth Tests of Goodness
%  of Fit", Oxford University Press, New York.
%
%  Scholz, F. W. and Stephens, M. A. (1987), "K-Sample
%  Anderson-Darling Tests", Journal of the American Statistical
%  Association, Vol. 82, No. 399, pp. 918-924.
%  DOI: 10.1080/01621459.1987.10478517.
%
%  Van Buuren, S. and Groothuis-Oudshoorn, K. (2011), "mice: Multivariate
%  Imputation by Chained Equations in R", Journal of Statistical Software,
%  Vol. 45, No. 3, pp. 1-67. DOI: 10.18637/jss.v045.i03.
%
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdJJtest')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
   %% Example 1: Jamshidian-Jalal test under MCAR.
   rng(10)
   n = 500;
   p = 5;
   Y = randn(n,p);
   Y(rand(n,p)<0.20) = NaN;
   out = mdJJtest(Y);
   disp(out.pvalue)
%}
%
%{
   %% Example 2: Display diagnostic plots.
   rng(20)
   n = 600;
   p = 6;
   Y = randn(n,p);
   Y(rand(n,p)<0.15) = NaN;
   out = mdJJtest(Y,'plots',true);
%}
%
%{
   %% Example 3: Use the nonparametric branch directly.
   rng(30)
   n = 500;
   p = 5;
   Y = trnd(5,n,p);
   Y(rand(n,p)<0.20) = NaN;
   out = mdJJtest(Y,'method','nonparametric');
   disp(out.ADpvalue)
%}
%
%{
   %% Example 4: Supply completed data sets.
   rng(40)
   n = 300;
   p = 4;
   Y = randn(n,p);
   Y(rand(n,p)<0.10) = NaN;
   outEM = mdEM(Y);
   Yimp = cell(3,1);
   for j = 1:3
       Yimp{j} = mdImputeStochastic(Y,outEM.loc,outEM.cov);
   end
   out = mdJJtest(Y,'imputed',Yimp,'method','hawkins');
%}

%{
    % Example 5: e input is a table.
    rng(50)
    Y = randn(300,4);
    Y(rand(size(Y))<0.15) = NaN;
    Ytable = array2table(Y,'VariableNames',{'Income','Age','Score','Balance'});
    out = mdJJtest(Ytable,'plots',true,'msg',false);
%}

%{
   %% Example 6: Hawkins rejection driven by one missingness pattern.
   % The first three missingness-pattern groups are generated from the same
   % multivariate standard normal distribution. The fourth and smallest
   % pattern has an inflated covariance matrix and is therefore expected to
   % dominate the Hawkins statistic and its Pareto plot.
   %
   % The rejection is not caused by the mere existence of the fourth
   % missingness pattern, but by the covariance heterogeneity associated
   % with that pattern. After removing it, the Hawkins test should no longer
   % be significant.
   rng(60)
   p = 5;
   groupSizes = [190; 190; 190; 30];
   group = repelem((1:4)',groupSizes);

   % Generate the complete data before introducing missing values.
   Ycomplete = randn(sum(groupSizes),p);
   Ycomplete(group==4,:) = 3*Ycomplete(group==4,:);

   % Introduce four distinct missingness patterns. Pattern P1 is the
   % complete-case pattern; P4 is the pattern with inflated covariance.
   Y = Ycomplete;
   Y(group==2,1) = NaN;
   Y(group==3,2) = NaN;
   Y(group==4,[3 4]) = NaN;

   % Supply the original complete data as the completed data set so that
   % the example focuses only on differences among missingness patterns.
   out = mdJJtest(Y,'imputed',Ycomplete,'method','hawkins', ...
       'plots',true,'msg',false);

   contributionTable = table((1:out.npatterns)', ...
       out.patternCounts,out.meanHawkinsContribution, ...
       'VariableNames',{'Pattern','Frequency','MeanContribution'});
   disp(contributionTable)
   fprintf('Hawkins p-value with all patterns: %.4g\n',out.pvalue)

   % Remove P4 and repeat the test. The remaining pattern groups have the
   % same distribution, so the rejection should disappear.
   keep = group~=4;
   outWithoutP4 = mdJJtest(Y(keep,:), ...
       'imputed',Ycomplete(keep,:),'method','hawkins', ...
       'plots',false,'msg',false);
   fprintf('Hawkins p-value without P4: %.4g\n',outWithoutP4.pvalue)
%}

%% Beginning of code

% The computational workflow is:
%
%   1. validate the incomplete data and the optional arguments;
%   2. identify the missingness patterns in the original data;
%   3. generate or validate one or more completed data sets;
%   4. remove patterns that are too sparse for reliable inference;
%   5. compute Hawkins' test on every completed data set;
%   6. when requested, compute the k-sample Anderson-Darling test;
%   7. summarize the imputation-specific results and construct the output.
%
% The missingness patterns are always determined from the original
% incomplete matrix Y. Imputed values are used only to evaluate the test
% statistics within those fixed pattern groups.

% Input parameters checking
if nargin < 1
    error('FSDA:mdJJtest:TooFewInputs', ...
        'At least one input argument is required.');
end

% Preserve table variable names for graphical output. The numerical
% calculations use only the table contents.
varNames = [];
if istable(Y)
    varNames = Y.Properties.VariableNames;
    Y = Y{:,:};
end

if ~ismatrix(Y) || ~isnumeric(Y)
    error('FSDA:mdJJtest:WrongInput', ...
        'Input argument Y must be a numeric matrix or a numeric table.');
end

% Work in double precision because covariance calculations, matrix
% factorizations and distribution functions are performed below.
Y = double(Y);
[nOriginal,p] = size(Y);

if nOriginal < 4 || p < 2
    error('FSDA:mdJJtest:TooSmall', ...
        'At least four observations and two variables are required.');
end

if any(all(isnan(Y),1))
    error('FSDA:mdJJtest:AllMissingColumn', ...
        'At least one variable contains only missing values.');
end

if any(all(isnan(Y),2))
    warning('FSDA:mdJJtest:AllMissingRow', ...
        ['Some rows contain only missing values. The test may be unreliable ' ...
        'for these observations.']);
end

if any(columnRangeIgnoringNaN(Y)==0)
    error('FSDA:mdJJtest:ConstantColumn', ...
        'At least one variable is constant on its observed values.');
end

% Default options. Options controlling imputation are used only when the
% completed data sets are not supplied explicitly.
imputed = [];
nimputations = 5;
minn = 6;
method = 'auto';
nsimul = 10000;
usechisq = 30;
alpha = 0.05;
maxiter = 200;
tol = 1e-7;
ridge = 1e-10;
msg = true;
plots = false;

% Parse optional name-value arguments. The structure options contains the
% canonical option names and their defaults. aux.chkoptions rejects unknown
% names before the user values are copied into the structure.
if ~isempty(varargin)
    options = struct('imputed',imputed,'nimputations',nimputations, ...
        'minn',minn,'method',method,'nsimul',nsimul, ...
        'usechisq',usechisq,'alpha',alpha, ...
        'maxiter',maxiter,'tol',tol,'ridge',ridge, ...
        'msg',msg,'plots',plots);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));

    if ~isempty(UserOptions)
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:mdJJtest:WrongInputOpt', ...
                ['Number of supplied options is invalid. Probably values ' ...
                'for some parameters are missing.']);
        end
        aux.chkoptions(options,UserOptions)
    end

    for j = 1:2:length(varargin)
        options.(varargin{j}) = varargin{j+1};
    end

    imputed = options.imputed;
    nimputations = options.nimputations;
    minn = options.minn;
    method = options.method;
    nsimul = options.nsimul;
    usechisq = options.usechisq;
    alpha = options.alpha;
    maxiter = options.maxiter;
    tol = options.tol;
    ridge = options.ridge;
    msg = options.msg;
    plots = options.plots;
end

% Basic checks on options. These checks are kept separate from parsing so
% that every option receives a method-specific and informative error.
if ~isempty(imputed) && ~(iscell(imputed) || isnumeric(imputed))
    error('FSDA:mdJJtest:WrongImputed', ...
        'Option ''imputed'' must be empty, a cell array or a numeric array.');
end

if ~isscalar(nimputations) || nimputations < 1 || ...
        nimputations ~= floor(nimputations)
    error('FSDA:mdJJtest:WrongNimputations', ...
        'Option ''nimputations'' must be a positive integer.');
end

if ~isnumeric(minn) || ~isscalar(minn) || ~isreal(minn) || ...
        ~isfinite(minn) || minn < 2 || minn ~= floor(minn)
    error('FSDA:mdJJtest:WrongMinn', ...
        ['Option ''minn'' must be a finite integer greater than or ' ...
        'equal to 2.']);
end

if ~(ischar(method) || (isstring(method) && isscalar(method)))
    error('FSDA:mdJJtest:WrongMethod', ...
        'Option ''method'' must be a character vector or a string scalar.');
end
method = lower(char(method));
if ~ismember(method,{'auto','hawkins','nonparametric'})
    error('FSDA:mdJJtest:WrongMethod', ...
        ['Option ''method'' must be ''auto'', ''hawkins'' or ' ...
        '''nonparametric''.']);
end

if ~isscalar(nsimul) || nsimul < 100 || nsimul ~= floor(nsimul)
    error('FSDA:mdJJtest:WrongNsimul', ...
        'Option ''nsimul'' must be an integer greater than or equal to 100.');
end

if ~isscalar(usechisq) || usechisq < 2 || usechisq ~= floor(usechisq)
    error('FSDA:mdJJtest:WrongUsechisq', ...
        'Option ''usechisq'' must be an integer greater than or equal to 2.');
end

if ~isscalar(alpha) || ~isnumeric(alpha) || alpha <= 0 || alpha >= 1
    error('FSDA:mdJJtest:WrongAlpha', ...
        'Option ''alpha'' must be a scalar in the interval (0,1).');
end

if ~isscalar(maxiter) || maxiter < 1 || maxiter ~= floor(maxiter)
    error('FSDA:mdJJtest:WrongMaxiter', ...
        'Option ''maxiter'' must be a positive integer.');
end

if ~isscalar(tol) || ~isnumeric(tol) || tol <= 0
    error('FSDA:mdJJtest:WrongTol', ...
        'Option ''tol'' must be a positive scalar.');
end

if ~isscalar(ridge) || ~isnumeric(ridge) || ridge < 0
    error('FSDA:mdJJtest:WrongRidge', ...
        'Option ''ridge'' must be a nonnegative scalar.');
end

if ~isscalar(msg) || ...
        ~(islogical(msg) || isnumeric(msg)) || ...
        ~(msg == 0 || msg == 1)
    error('FSDA:mdJJtest:WrongMsg', ...
        'Option ''msg'' must be true, false, 0, or 1.');
end

if ~isscalar(plots) || ...
        ~(islogical(plots) || isnumeric(plots)) || ...
        ~(plots == 0 || plots == 1)
    error('FSDA:mdJJtest:WrongPlots', ...
        'Option ''plots'' must be true, false, 0, or 1.');
end

msg = logical(msg);
plots = logical(plots);


% Identify the original missingness patterns. A true entry means that the
% corresponding variable is missing. Rows with identical logical patterns
% are assigned the same integer label in idxPatternsOriginal.
missingMaskOriginal = isnan(Y);
[patternsOriginal,idxPatternsOriginal,countsOriginal] = ...
    patternGroups(missingMaskOriginal);

% Construct or validate completed data sets. The same original pattern
% labels will be used for every completed data set, so the test compares
% groups defined by the observed missingness structure rather than by the
% imputed values.
loc = [];
covmat = [];
if isempty(imputed)
    % Estimate the joint-normal location and covariance using all
    % observations and the already computed pattern information.
    outEM = mdEM(Y,'maxiter',maxiter,'tol',tol, ...
        'Patterns',patternsOriginal,'idxPatterns',idxPatternsOriginal);
    loc = outEM.loc;
    covmat = outEM.cov;
    % Generate independent stochastic completions conditional on the EM
    % estimates. Each completion is analysed separately below.
    imputations = cell(nimputations,1);
    for j = 1:nimputations
        imputations{j} = mdImputeStochastic(Y,loc,covmat);
    end
else
    imputations = normalizeImputedInput(imputed,Y);
end

% Remove sparse patterns. A pattern is retained only if its frequency is
% greater than or equal to minn. The removal is based on the pattern
% frequencies in the original incomplete data and is applied identically to
% every completed data set.
removePattern = countsOriginal < minn;
removedRows = removePattern(idxPatternsOriginal);
removedPatterns = patternsOriginal(removePattern,:);
% removedPatternCounts = countsOriginal(removePattern);

if all(removedRows)
    error('FSDA:mdJJtest:NoRows', ...
        ['All observations belong to missingness patterns with frequency ' ...
        'less than minn. Lower minn with caution.']);
end

% Restrict both the incomplete data and all completed data sets to the same
% retained observations.
Yused = Y(~removedRows,:);
imputations = cellfun(@(Z) Z(~removedRows,:),imputations, ...
    'UniformOutput',false);

% Recompute the pattern encoding after sparse rows have been removed.
% patterns(r,:) describes retained group r and idxPatterns(i) gives the
% group membership of retained observation i.
[patterns,idxPatterns,patternCounts] = patternGroups(isnan(Yused));
g = size(patterns,1);
nused = size(Yused,1);

if g < 2
    error('FSDA:mdJJtest:OnePattern', ...
        'At least two retained missingness patterns are required.');
end

if any(patternCounts < 2)
    error('FSDA:mdJJtest:TinyPattern', ...
        'At least two observations are required in every retained pattern.');
end

if nused-g-p <= 0
    error('FSDA:mdJJtest:InsufficientDF', ...
        'Hawkins'' transformation requires nused-g-p to be positive.');
end

% Compute Hawkins' transformed quantities for every completed data set.
% details{j} stores the pattern-specific transformed values for imputation j
% and is also reused if the nonparametric branch is required.
M = numel(imputations);
details = cell(M,1);
for j = 1:M
    details{j} = hawkinsCore(imputations{j},idxPatterns,ridge);
end

hawkinsStat = [];
hawkinsP = [];
hawkinsDF = 2*g;

% Hawkins' test is required for methods 'auto' and 'hawkins'. For every
% imputation, one uniformity p-value is computed for each retained pattern
% and the pattern-specific p-values are combined using Fisher's method.
if strcmp(method,'auto') || strcmp(method,'hawkins')
    hawkinsStat = zeros(M,1);
    hawkinsP = zeros(M,1);

    for j = 1:M
        groupP = zeros(g,1);
        for r = 1:g
            groupP(r) = neymanPValue(details{j}.A{r},nsimul,usechisq);
            if groupP(r)==0
                groupP(r) = 1/nsimul;
            end
        end

        % Fisher's combination statistic is asymptotically chi-square with
        % 2*g degrees of freedom when the pattern-specific p-values are
        % independent under the null.
        hawkinsStat(j) = -2*sum(log(groupP));
        hawkinsP(j) = chi2cdf(hawkinsStat(j),hawkinsDF,"upper");
        details{j}.NeymanP = groupP;
    end
end

meanHawkinsContribution = [];

if ~isempty(hawkinsP)
    meanHawkinsContribution = zeros(g,1);

    for j = 1:M
        hawkinsContribution = -2*log(details{j}.NeymanP(:));
        details{j}.HawkinsContributions = hawkinsContribution;
        meanHawkinsContribution = meanHawkinsContribution + ...
            hawkinsContribution;
    end

    meanHawkinsContribution = meanHawkinsContribution/M;
end

% Compute the k-sample Anderson-Darling test when required. In automatic
% mode this branch is entered when at least one imputation-specific Hawkins
% p-value is below alpha, following the logic used in mice::mcar.
ADstat = [];
ADpvalue = [];
meanADContribution = [];

runAD = strcmp(method,'nonparametric') || ...
    (strcmp(method,'auto') && any(hawkinsP<alpha));

if runAD
    ADstat = zeros(M,1);
    ADpvalue = zeros(M,1);
    meanADContribution = zeros(g,1);

    for j = 1:M
        ADout = andersonDarlingCore(details{j}.Fij);

        ADstat(j) = ADout.Statistic;
        ADpvalue(j) = ADout.PValue;

        details{j}.AndersonDarling = ADout;
        meanADContribution = meanADContribution + ...
            ADout.GroupContributions(:);
    end

    meanADContribution = meanADContribution/M;
end



% Select the statistic reported in out.stat and out.pvalue. Multiple
% imputations are summarized by their medians, while the full vectors remain
% available in out.hawkinsP, out.ADpvalue and the corresponding statistics.
if strcmp(method,'hawkins') || (strcmp(method,'auto') && isempty(ADpvalue))
    stat = median(hawkinsStat);
    pvalue = median(hawkinsP);
    df = hawkinsDF;
else
    stat = median(ADstat);
    pvalue = median(ADpvalue);
    df = [];
end

% Create the pattern-information table. Each row name is the binary
% missingness pattern written in variable order, where 1 denotes missing and
% 0 denotes observed.
patternNames = cell(g,1);
for r = 1:g
    patternNames{r} = sprintf('%d',patterns(r,:));
end
patternInfo = table(patternCounts,sum(patterns,2), ...
    'VariableNames',{'Frequency','NMissing'}, ...
    'RowNames',patternNames);

% Store results. Both the aggregate decision quantities and the
% imputation-specific diagnostics are returned so that the user can inspect
% sensitivity to the stochastic completions.
out = struct;
out.stat = stat;
out.pvalue = pvalue;
out.df = df;
out.method = method;
out.alpha = alpha;
out.minn = minn;
out.nsimul = nsimul;
out.usechisq = usechisq;
out.nimputations = M;
out.n = nOriginal;
out.nused = nused;
out.p = p;
out.patterns = patterns;
out.patternCounts = patternCounts;
out.patternInfo = patternInfo;
out.npatterns = g;
out.idxPatterns = idxPatterns;
out.removedRows = removedRows;
out.removedPatterns = removedPatterns;
out.hawkinsStat = hawkinsStat;
out.hawkinsDF = hawkinsDF;
out.hawkinsP = hawkinsP;
out.ADstat = ADstat;
out.ADpvalue = ADpvalue;
out.loc = loc;
out.cov = covmat;
out.details = details;
out.interpretation = interpretResult(method,hawkinsP,ADpvalue,alpha);


out.meanHawkinsContribution = meanHawkinsContribution;
out.meanADContribution = meanADContribution;

if msg
    printSummary(out)
end

if plots
    plotResults(out,Y,varNames)
end
end

% -------------------------------------------------------------------------
function imputations = normalizeImputedInput(imputed,Y)
%normalizeImputedInput Convert supplied completed data sets to a cell array.
%
% The function accepts one completed matrix, an n-by-p-by-M numeric array,
% or a cell array of completed matrices/tables. It also verifies that the
% observed entries of Y have not been altered, apart from numerical roundoff.

[n,p] = size(Y);

if iscell(imputed)
    imputations = imputed(:);
elseif isnumeric(imputed) && ndims(imputed)<=3
    if ismatrix(imputed)
        imputations = {double(imputed)};
    else
        M = size(imputed,3);
        imputations = cell(M,1);
        for j = 1:M
            imputations{j} = double(imputed(:,:,j));
        end
    end
else
    error('FSDA:mdJJtest:WrongImputed', ...
        ['Option ''imputed'' must be a cell array or an n x p x M ' ...
        'numeric array.']);
end

% The comparison with the original data is restricted to entries that were
% observed before imputation. Missing entries may of course differ across
% completed data sets.
observed = ~isnan(Y);
observedValues = Y(observed);
if isempty(observedValues)
    toleranceObserved = 0;
else
    toleranceObserved = 100*eps(max(1,max(abs(observedValues))));
end

% Validate every completed data set separately.
for j = 1:numel(imputations)
    if istable(imputations{j})
        imputations{j} = imputations{j}{:,:};
    end

    if ~isnumeric(imputations{j}) || ...
            ~isequal(size(imputations{j}),[n p])
        error('FSDA:mdJJtest:WrongImputedSize', ...
            'Every completed data set must have the same size as Y.');
    end

    imputations{j} = double(imputations{j});

    if any(isnan(imputations{j}(:)))
        error('FSDA:mdJJtest:IncompleteImputation', ...
            'A supplied completed data set still contains NaN values.');
    end

    if any(abs(imputations{j}(observed)-Y(observed))>toleranceObserved)
        warning('FSDA:mdJJtest:ObservedChanged', ...
            'Supplied imputation %d changes at least one observed value.',j);
    end
end
end

% -------------------------------------------------------------------------
function [patterns,idxPatterns,counts] = patternGroups(missingMask)
%patternGroups Group rows having the same missingness pattern.
%
% patterns contains the distinct rows in order of first appearance;
% idxPatterns maps each observation to a row of patterns; counts contains the
% corresponding group frequencies.

[patterns,~,idxPatterns] = unique(missingMask,'rows','stable');
counts = accumarray(idxPatterns,1,[size(patterns,1) 1]);
end

% -------------------------------------------------------------------------
function out = hawkinsCore(Y,idxPatterns,ridge)
%hawkinsCore Compute Hawkins' transformed quantities.
%
% Y is one completed data set. Observations are partitioned according to
% idxPatterns, which was obtained from the original incomplete data. The
% function first estimates a covariance matrix pooled across pattern groups
% and then transforms each within-group squared Mahalanobis distance to an F
% quantity and, finally, to a value that should be Uniform(0,1) under the
% Hawkins null model.

[n,p] = size(Y);
groups = unique(idxPatterns,'stable');
g = length(groups);

% pooledSS accumulates the within-pattern sums of squares and cross-products.
% Group means are removed separately so that between-pattern mean differences
% do not contribute to the pooled covariance estimate.
pooledSS = zeros(p,p);
centeredGroups = cell(g,1);
groupSizes = zeros(g,1);

for r = 1:g
    Yr = Y(idxPatterns==groups(r),:);
    groupSizes(r) = size(Yr,1);
    Yr = Yr-mean(Yr,1);
    centeredGroups{r} = Yr;
    pooledSS = pooledSS+Yr'*Yr;
end

% There are g separately estimated group means, giving n-g residual degrees
% of freedom for the pooled covariance matrix.
pooledCov = pooledSS/(n-g);
pooledCov = stabilizeSPD(pooledCov,ridge,'pooled covariance');
invPooledCov = pooledCov\eye(p);
invPooledCov = (invPooledCov+invPooledCov')/2;

Fij = cell(g,1);
A = cell(g,1);
df2 = n-g-p;

for r = 1:g
    Yr = centeredGroups{r};

    % Row-wise squared Mahalanobis distances based on the pooled covariance.
    mahal = sum((Yr*invPooledCov).*Yr,2);
    h = groupSizes(r)*mahal;
    den = p*((groupSizes(r)-1)*(n-g)-h);

    if any(den<=0)
        error('FSDA:mdJJtest:HawkinsDenominator', ...
            ['A nonpositive denominator occurred in Hawkins'' ' ...
            'transformation. Check pattern sizes, collinearity or the ' ...
            'supplied imputations.']);
    end

    % Hawkins' finite-sample transformation. Under multivariate normality
    % and equality of covariance matrices across patterns, Fij follows the
    % stated F reference law and its upper-tail probability A should be
    % approximately Uniform(0,1).
    Fij{r} = ((n-g-p)*h)./den;
    A{r} = fSurvival(Fij{r},p,df2);
end

out = struct;
out.Fij = Fij;
out.A = A;
out.groupSizes = groupSizes;
out.pooledCov = pooledCov;
end

% -------------------------------------------------------------------------
function pvalue = neymanPValue(u,nsimul,usechisq)
%neymanPValue Fourth-order Neyman smooth test for uniformity.
%
% The test projects the empirical distribution of u onto the first four
% orthonormal shifted Legendre polynomials. Under exact uniformity the four
% normalized projection coefficients are asymptotically standard normal, so
% their squared sum has a chi-square distribution with four degrees of
% freedom. Small samples are calibrated by Monte Carlo simulation instead.

u = u(:);
n = length(u);
% Phi(i,k) is the kth orthonormal polynomial evaluated at u(i). The column
% sums measure departures from uniformity in four increasingly complex
% directions.
Phi = shiftedLegendre(u,4);
stat = sum(sum(Phi,1).^2)/n;

if n<usechisq
    % For small pattern sizes, simulate the null distribution using nsimul
    % independent samples of size n from Uniform(0,1).
    U = rand(n,nsimul);
    z = 2*U-1;
    P1 = z;
    P2 = (3*z.^2-1)/2;
    P3 = (5*z.*P2-2*P1)/3;
    P4 = (7*z.*P3-3*P2)/4;
    sumPol = [sqrt(3)*sum(P1,1); sqrt(5)*sum(P2,1); ...
        sqrt(7)*sum(P3,1); 3*sum(P4,1)];
    simulated = sum(sumPol.^2,1)/n;
    pvalue = mean(simulated>stat);
else
    % For sufficiently large patterns, use the asymptotic chi-square law.
    pvalue = chi2cdf(stat,4,"upper");
end
end

% -------------------------------------------------------------------------
function Phi = shiftedLegendre(u,order)
%shiftedLegendre Orthonormal shifted Legendre polynomials.
%
% Standard Legendre polynomials are evaluated at z=2*u-1, which maps the
% interval [0,1] to [-1,1]. Multiplication by sqrt(2*degree+1) makes the
% resulting shifted polynomials orthonormal under the Uniform(0,1) measure.

u = u(:);
z = 2*u-1;
Phi = zeros(length(u),order);
Pprev = ones(size(z));
Pcurr = z;
Phi(:,1) = sqrt(3)*Pcurr;

% Three-term Legendre recurrence, avoiding repeated calls to symbolic or
% polynomial-construction routines.
for degree = 2:order
    Pnext = ((2*degree-1)*z.*Pcurr-(degree-1)*Pprev)/degree;
    Phi(:,degree) = sqrt(2*degree+1)*Pnext;
    Pprev = Pcurr;
    Pcurr = Pnext;
end
end

% -------------------------------------------------------------------------
function out = andersonDarlingCore(Fij)
%andersonDarlingCore k-sample Anderson-Darling test.
%
% Fij is a cell array containing one sample for each retained missingness
% pattern. The test compares the complete empirical distributions across
% patterns and is therefore robust to common nonnormal shapes that may make
% Hawkins' normal-theory test significant. The implementation follows the
% Scholz-Stephens k-sample Anderson-Darling approximation used by mice.

allF = vertcat(Fij{:});
groupSizes = cellfun(@numel,Fij);
k = length(groupSizes);
n = length(allF);

if k<2
    error('FSDA:mdJJtest:ADGroups', ...
        'At least two groups are required for the Anderson-Darling test.');
end

if n<4
    error('FSDA:mdJJtest:ADSample', ...
        'At least four observations are required for the Anderson-Darling test.');
end

% The largest pooled observation is excluded because the Anderson-Darling
% denominator contains H_j*(n-H_j), which is zero at the final order statistic.
xsort = sort(allF);
xsort(end) = [];
isNew = [true; diff(xsort)~=0];
zj = xsort(isNew);
first = find(isNew);
last = [first(2:end)-1; length(xsort)];
hj = last-first+1;
hn = cumsum(hj);

% Accumulate the contribution of each group over the distinct pooled order
% statistics. Ties are handled through hj, the multiplicity of each value.
ADgroupRaw = zeros(k,1);
for r = 1:k
    fr = Fij{r}(:);
    counts = zeros(length(zj),1);
    for j = 1:length(zj)
        counts(j) = sum(fr==zj(j));
    end
    mij = cumsum(counts);
    num = (n*mij-groupSizes(r)*hn).^2;
    den = hn.*(n-hn);
    ADgroupRaw(r) = sum(hj.*(num./den))/groupSizes(r);
end

% Overall k-sample Anderson-Darling statistic and normalized group
% contributions.
ADstat = sum(ADgroupRaw)/n;
ADgroup = ADgroupRaw/n;
J = sum(1./groupSizes);
harmonic = cumsum(1./(1:(n-1)));
H = harmonic(end);
idx = 1:(n-2);
G = sum((H-harmonic(idx))./(n-idx));

% Finite-sample variance approximation for the standardized statistic.
% H and G are harmonic-sum quantities appearing in the k-sample theory.
a = (4*G-6)*(k-1)+(10-6*G)*J;
b = (2*G-4)*k^2+8*H*k+(2*G-14*H-4)*J-8*H+4*G-6;
c = (6*H+2*G-2)*k^2+(4*H-4*G+6)*k+(2*H-6)*J+4*H;
d = (2*H+6)*k^2-4*H*k;
variance = max(((a*n^3)+(b*n^2)+(c*n)+d)/ ...
    ((n-1)*(n-2)*(n-3)),0);

if variance==0
    standardized = sign(ADstat-(k-1))*Inf;
else
    standardized = (ADstat-(k-1))/sqrt(variance);
end

% Approximate the upper-tail probability by interpolating tabulated
% standardized quantiles on the logit scale.
b0 = [0.675 1.281 1.645 1.960 2.326];
b1 = [-0.245 0.250 0.678 1.149 1.822];
b2 = [-0.105 -0.305 -0.362 -0.391 -0.396];
c0 = [1.09861228866811 2.19722457733622 2.94443897916644 ...
    3.66356164612965 4.59511985013459];
qnt = b0+b1/sqrt(k-1)+b2/(k-1);

if standardized<=qnt(3)
    take = 2:5;
else
    take = 1:4;
end

yy = spline(qnt(take),c0(take),standardized);
if yy>=0
    ey = exp(-yy);
    pvalue = ey/(1+ey);
else
    pvalue = 1/(1+exp(yy));
end

out = struct;
out.PValue = pvalue;
out.Statistic = ADstat;
out.GroupContributions = ADgroup;
out.Standardized = standardized;
out.Variance = variance;
end

% -------------------------------------------------------------------------
function A = stabilizeSPD(A,ridge,label)
%stabilizeSPD Symmetrize and regularize a covariance matrix if necessary.
%
% A small diagonal ridge is always added when ridge is positive. If Cholesky
% factorization still fails, the eigenvalues are floored at a scale-dependent
% positive value, producing the nearest matrix used here on the same
% eigenvector basis.

A = (A+A')/2;
if isempty(A)
    return
end

base = max(1,trace(abs(A))/size(A,1));
if ridge>0
    A = A+ridge*base*eye(size(A));
end

[~,flag] = chol(A);
if flag~=0
    [V,D] = eig(A);
    d = real(diag(D));
    floorValue = max(ridge*base,eps(base));
    d = max(d,floorValue);
    A = V*diag(d)*V';
    A = (A+A')/2;
    warning('FSDA:mdJJtest:Regularized', ...
        'The %s was projected onto the positive-definite cone.',label);
end
end

% -------------------------------------------------------------------------
function pvalue = fSurvival(x,df1,df2)
%fSurvival Upper-tail probability of the F distribution.
%
% The beta-function representation avoids requiring a direct call to fcdf
% and is numerically stable for the nonnegative statistics used here.

x = max(x,0);
z = df2./(df2+df1*x);
pvalue = betainc(z,df2/2,df1/2);
pvalue(isinf(x)) = 0;
end


% -------------------------------------------------------------------------
function txt = interpretResult(method,hawkinsP,ADpvalue,alpha)
%interpretResult Produce a concise interpretation of the selected test.
%
% The wording distinguishes three situations: no evidence against the joint
% Hawkins null, evidence against MCAR from the nonparametric branch, and a
% Hawkins rejection that is plausibly attributable to nonnormality.

if strcmp(method,'hawkins')
    if median(hawkinsP)<alpha
        txt = ['Hawkins'' test is significant. Under multivariate ' ...
            'normality this is evidence against MCAR.'];
    else
        txt = ['Hawkins'' test is not significant: there is no evidence ' ...
            'against multivariate normality or MCAR.'];
    end
elseif strcmp(method,'nonparametric')
    if median(ADpvalue)<alpha
        txt = ['The Anderson-Darling rank test is significant: reject ' ...
            'MCAR within the MCAR-versus-MAR framework.'];
    else
        txt = ['The Anderson-Darling rank test is not significant: ' ...
            'there is no evidence against MCAR.'];
    end
else
    if isempty(ADpvalue)
        txt = ['No Hawkins p-value is smaller than alpha: there is no ' ...
            'evidence against multivariate normality or MCAR.'];

    elseif median(ADpvalue)<alpha
        txt = ['At least one Hawkins p-value is smaller than alpha and ' ...
            'the median Anderson-Darling p-value is significant. This ' ...
            'provides evidence against MCAR within the MCAR-versus-MAR ' ...
            'framework.'];

    else
        txt = ['At least one Hawkins p-value is smaller than alpha, but ' ...
            'the median Anderson-Darling p-value is not significant. ' ...
            'This suggests nonnormality, with no evidence against MCAR.'];
    end
end
end

% -------------------------------------------------------------------------
function printSummary(out)
%printSummary Display a concise summary in the Command Window.

fprintf('\nJamshidian-Jalal MCAR test\n')
fprintf('Retained patterns: %d\n',size(out.patterns,1))
fprintf('Cases used: %d of %d\n',out.nused,out.n)

if any(out.removedRows)
    fprintf('Cases removed because of sparse patterns: %d\n', ...
        sum(out.removedRows))
end

if ~isempty(out.hawkinsP)
    fprintf(['Hawkins: median chi-square = %.6g, df = %d, ' ...
        'median p-value = %.6g\n'],median(out.hawkinsStat), ...
        out.hawkinsDF,median(out.hawkinsP))

    if any(out.hawkinsP<out.alpha) && median(out.hawkinsP)>=out.alpha
        fprintf('Warning: at least one Hawkins p-value is below alpha.\n')
    end
end

if ~isempty(out.ADpvalue)
    fprintf(['Anderson-Darling: median statistic = %.6g, ' ...
        'median p-value = %.6g\n'],median(out.ADstat), ...
        median(out.ADpvalue))

    if any(out.ADpvalue<out.alpha) && median(out.ADpvalue)>=out.alpha
        fprintf(['Warning: at least one Anderson-Darling p-value is ' ...
            'below alpha.\n'])
    end
end

fprintf('%s\n',out.interpretation)
fprintf(['Caution: the procedure does not distinguish MCAR from MNAR; ' ...
    'a nonsignificant result does not rule out MNAR.\n\n'])
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function plotResults(out,Y,varNames)
%plotResults Display p-values, retained patterns and contributions.
%
% The first figure summarizes the variation of p-values across completed
% data sets. The second figure displays the retained missingness patterns.
% The third figure contains Pareto plots of the mean pattern contributions
% to the available test statistics.

ntests = (~isempty(out.hawkinsP)) + (~isempty(out.ADpvalue));

%% P-values across completed data sets
if ntests>0
    figure('Name','mdJJtest p-values','Color','w')
    tl = tiledlayout(ntests,1, ...
        'TileSpacing','compact','Padding','compact');

    if ~isempty(out.hawkinsP)
        ax = nexttile(tl);

        medianP = median(out.hawkinsP);

        bar(ax,out.hawkinsP)
        yline(ax,out.alpha,'--')
        xlabel(ax,'Completed data set')
        ylabel(ax,'Hawkins p-value')

        title(ax,sprintf(['Median p-value = %.4g; ' ...
            '%.1f%% below alpha'], ...
            medianP,100*mean(out.hawkinsP<out.alpha)))
    end

    if ~isempty(out.ADpvalue)
        ax = nexttile(tl);

        medianP = median(out.ADpvalue);

        bar(ax,out.ADpvalue)
        yline(ax,out.alpha,'--')
        xlabel(ax,'Completed data set')
        ylabel(ax,'Anderson-Darling p-value')

        title(ax,sprintf(['Median p-value = %.4g; ' ...
            '%.1f%% below alpha'], ...
            medianP,100*mean(out.ADpvalue<out.alpha)))
    end
end

%% Missingness patterns retained for the MCAR test
Yretained = Y(~out.removedRows,:);
% Let mdpattern choose automatically between the balloon plot and the
% heatmap. Suppress the explanatory text normally displayed in the
% Command Window.
plo = struct;
plo.showExplanation = false;

% Construct labels in the same pattern order used internally by mdpattern.
% The permanent identifiers P1, P2, ... refer to the rows of out.patterns
% and are therefore the same identifiers used in the Pareto plots.
missingMask = isnan(Yretained);
[~,ordercols] = sort(sum(missingMask,1),'ascend');
patternsMd = ~flipud(unique(~missingMask(:,ordercols),'rows'));
patternsMdOriginalOrder = false(size(patternsMd));
patternsMdOriginalOrder(:,ordercols) = patternsMd;

[matched,patternNumber] = ismember(patternsMdOriginalOrder, ...
    out.patterns,'rows');
if ~all(matched)
    error('FSDA:mdJJtest:PatternLabelMismatch', ...
        'Unable to match the mdpattern rows to the retained patterns.');
end

plo.rowLabels = compose('P%d (n=%d)',patternNumber, ...
    out.patternCounts(patternNumber));

if isempty(varNames)
    mdpattern(Yretained,'plots',plo);
else
    mdpattern(Yretained,'Lc',varNames,'plots',plo);
end

set(gcf,'Name','Retained missingness patterns')

%% Mean pattern contributions
ncontributions = ...
    (~isempty(out.meanHawkinsContribution)) + ...
    (~isempty(out.meanADContribution));

if ncontributions>0
    figure('Name','Pattern contributions','Color','w')
    tl = tiledlayout(ncontributions,1, ...
        'TileSpacing','compact','Padding','compact');

    if ~isempty(out.meanHawkinsContribution)
        ax = nexttile(tl);
        contributionPareto(ax,out.meanHawkinsContribution, ...
            out.patternCounts,'Hawkins');
    end

    if ~isempty(out.meanADContribution)
        ax = nexttile(tl);
        contributionPareto(ax,out.meanADContribution, ...
            out.patternCounts,'Anderson-Darling');
    end
end

end


% -------------------------------------------------------------------------
function ranges = columnRangeIgnoringNaN(Y)
%columnRangeIgnoringNaN Range of each variable over its observed values.

ranges = zeros(1,size(Y,2));
for j = 1:size(Y,2)
    values = Y(~isnan(Y(:,j)),j);
    ranges(j) = max(values)-min(values);
end
end


% -------------------------------------------------------------------------
function contributionPareto(ax,meanContribution,patternCounts,testName)
%contributionPareto Pareto plot of mean contributions by pattern.
%
% Bars show mean additive contributions in decreasing order. The line on
% the right axis reports the cumulative percentage of the total mean
% contribution.

meanContribution = meanContribution(:);
patternCounts = patternCounts(:);
g = numel(meanContribution);

% Sort patterns from largest to smallest contribution.
[sortedContribution,ord] = sort(meanContribution,'descend');

% Cumulative percentage of the mean statistic.
totalContribution = sum(sortedContribution);

if totalContribution>0
    cumulativePercentage = ...
        100*cumsum(sortedContribution)/totalContribution;
else
    cumulativePercentage = zeros(g,1);
end

% Labels preserve the original pattern numbers.
patternLabels = compose('P%d (n=%d)', ...
    ord(:),patternCounts(ord));

yyaxis(ax,'left')
bar(ax,1:g,sortedContribution)
ylabel(ax,'Mean contribution')

yyaxis(ax,'right')
plot(ax,1:g,cumulativePercentage,'-o', ...
    'LineWidth',1.25,'MarkerSize',4)
yline(ax,80,'--','80%')
ylabel(ax,'Cumulative contribution (%)')
ylim(ax,[0 105])
yticks(ax,0:20:100)

xlim(ax,[0.5 g+0.5])
xticks(ax,1:g)
xticklabels(ax,patternLabels)
xtickangle(ax,45)
xlabel(ax,'Missingness pattern')

title(ax,sprintf(['Mean contribution to the %s statistic ' ...
    'across completed data sets'],testName))

grid(ax,'on')

end
%FScategory:MULT-MissingData
