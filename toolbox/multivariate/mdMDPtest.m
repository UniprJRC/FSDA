function out = mdMDPtest(Y, varargin)
%mdMDPtest Mahalanobis Distance Perturbation test for MCAR
%
%<a href="matlab: docsearchFS('mdMDPtest')">Link to the help function</a>
%
%  mdMDPtest implements the Mahalanobis Distance Perturbation (MDP) test
%  for assessing the Missing Completely At Random (MCAR) assumption. The
%  test is based on the change in Mahalanobis distances for the units
%  without missing values when location and scatter are estimated:
%
%    1) using only the complete rows;
%    2) using all rows through EM/TEM in the presence of missing values.
%
%  Mean-statistic inference is based on the analytical asymptotic reference
%  whenever the corresponding analytical branch is defined. If bootstrap=true,
%  an additional mask-preserving full-pipeline bootstrap calibration is
%  computed. The same bootstrap samples also calibrate the Hausman-type omnibus
%  statistic whenever its scatter-discrepancy covariance is available; the
%  discrepancy sandwich and statistical spectral rank are recomputed inside
%  every retained replicate. The null generator is selected through
%  bootstraptype: the backward-compatible default is a Gaussian model fitted
%  on the complete rows, while 'empiricalradial' uses a null-enforcing
%  empirical-radial elliptical generator. In both cases the observed
%  missingness mask is imposed on each generated sample and the complete MDP
%  pipeline is rerun.
%
%
%  Required input arguments:
%
%    Y : Input data. Matrix. n x p data matrix possibly containing NaNs.
%        Rows of Y represent observations and columns represent variables.
%        Data Types - double
%
%
%  Optional input arguments:
%
%    alpha : Trimming/robustness level. Scalar in the interval [0,0.5].
%            The default value is 0.
%
%            If alpha=0, classical complete-case estimates are used and
%            mdEM is used for the incomplete-data fit.
%
%            If alpha>0 and robust='FS', the Forward Search on the
%            complete cases starts monitoring at
%
%                h = floor(nComplete*(1-alpha)).
%
%            where nComplete is the number of rows with complete
%            observations.
%
%            The final complete-case location and scatter are not based
%            necessarily on exactly h observations. FSM applies its
%            signal-detection and stopping rules and the final estimates
%            are computed from the observations not declared as outliers.
%
%            If alpha>0 and robust='MCD', alpha is used as the MCD
%            breakdown-point option. If robust='MM', alpha is used as the
%            breakdown point of the preliminary S estimator.
%
%            For all robust estimators, outliers are declared using a
%            Bonferronized simultaneous confidence level of 0.99. For FSM
%            two signal/stopping rules are available. With bonflev=0.99,
%            the direct 99 percent Bonferroni bound is used. With
%            bonflev=[], the standard FSM signal-detection, validation and
%            envelope-resuperimposition rules are used. For MCD and MM,
%            the individual confidence level is 1-0.01/nComplete.
%            Example - 'alpha', 0.1
%            Data Types - double inside [0 0.5]
%
%
%   method : Rescaling method used inside mdTEM and mdPartialMD2full.
%            The default value is 'pri'. The analytical alpha>0 sandwich is
%            available for 'pri', 'expScale', 'zMap', 'chiMap' and
%            'betaMap'. The option 'detMap' is deliberately excluded from
%            the analytical sandwich because its distance adjustment depends
%            explicitly on the scatter
%            matrix and would require additional derivative terms.
%            Example - 'method','betaMap'
%            Data Types - char | string
%
% consistencyfactor : TEM truncation-bias correction used when alpha>0.
%            Character vector or string scalar. Possible values are
%            'adaptive', 'pattern', 'global', 'weighted' or 'none'. The
%            default is 'adaptive'.
%
%            With 'adaptive', mdTEM estimates radial consistency factors
%            from the complete cases, using the complete-case location and
%            scatter fitted by the selected robust estimator. Thus the
%            all-data TEM scatter is calibrated to the same affine-
%            equivariant elliptical scatter target as the complete-case
%            estimator without specifying a parametric radial distribution.
%            In particular, a Student t radial law does not require the
%            degrees of freedom to be known or estimated for this matching.
%            The reference sample and reference fit are reconstructed inside
%            every bootstrap replicate.
%
%            'pattern' uses Gaussian pattern-wise Tallis factors and is
%            retained as the Gaussian-reference correction. The options
%            'global', 'weighted' and 'none' have the same meaning as in
%            mdTEM. This option is ignored when alpha=0. The current
%            coupledtrim sensitivity implementation supports only
%            consistencyfactor='pattern'.
%            Example - 'consistencyfactor','adaptive'
%            Data Types - char | string
%
% adaptivepool : Pool adaptive projected reference radii across patterns
%            having the same observed dimension. Logical scalar. The
%            default is true. This option is used only when alpha>0 and
%            consistencyfactor='adaptive'. Under ellipticity, patterns of
%            the same observed dimension have the same standardized radial
%            law. Pooling stabilizes the factor estimate while the analytical
%            sandwich retains the dependence among the several projections
%            generated by each complete observation.
%            Example - 'adaptivepool',false
%            Data Types - logical
%
% adaptiveminref : Minimum number of complete-case reference radii inside
%            a pattern cutoff required by the adaptive correction. Positive
%            integer. The default is 20. This option is used only when
%            alpha>0 and consistencyfactor='adaptive'.
%            Example - 'adaptiveminref',10
%            Data Types - double
%
% omnibusrankexp : Exponent used by the statistical spectral rank rule for
%            the Hausman-type omnibus statistic. Scalar in the open interval
%            (0,0.5). The default is 1/3. If lambda_1 is the largest estimated
%            eigenvalue of the scatter-discrepancy covariance, the relative
%            threshold is
%
%                eta_n = nEff^(-omnibusrankexp)
%
%            and eigenvalue j is retained when lambda_j/lambda_1 > eta_n.
%            The values 1/4, 1/3 and 2/5 are also reported as sensitivity
%            diagnostics in out.omnibus.rankSensitivity.
%            Example - 'omnibusrankexp',0.4
%            Data Types - double
%
% omnibusneff : Effective sample size used by the statistical spectral rank
%            threshold. Character vector, string scalar or numeric scalar.
%            The default is 'n'. Possible character values are 'n' and
%            'complete'. With 'n', nEff is the total number of observations.
%            With 'complete', nEff is the number of complete cases. A numeric
%            scalar greater than one can be supplied directly. The option
%            affects only omnibus spectral rank selection; it does not change
%            fitting, scalar MDP statistics or the discrepancy covariance.
%            Example - 'omnibusneff','complete'
%            Data Types - char | string | double
%
% bootstrap : Add mask-preserving full-pipeline bootstrap calibration.
%            Logical scalar. The default is false. The analytical p-values
%            for the two mean statistics are evaluated independently of this
%            option whenever the corresponding analytical branch is defined.
%            If bootstrap=true, bootstrap p-values for all four scalar
%            statistics and bootstrap confidence intervals are additionally
%            computed. When the Hausman-type omnibus covariance is available,
%            a full-pipeline bootstrap p-value is also computed. In every
%            replicate the complete-case fit, the all-data EM/TEM fit, the
%            adaptive radial correction when requested, the selected mean
%            aggregation rule, the scatter-discrepancy covariance and the
%            omnibus spectral rank are recomputed.
%            Example - 'bootstrap',true
%            Data Types - logical
%
% bootstraptype : Bootstrap null generator. Character vector or string scalar.
%            Possible values are 'gaussian' and 'empiricalradial'. The
%            default is 'gaussian' for backward compatibility. This option
%            is used only when bootstrap=true.
%
%            'gaussian' generates full p-dimensional observations from the
%            Gaussian model fitted to the complete cases and then imposes the
%            observed missingness mask.
%
%            'empiricalradial' is a null-enforcing semiparametric elliptical
%            generator. Complete-case squared Mahalanobis radii are computed
%            in the fitted complete-case geometry and resampled with
%            replacement. The empirical radial law uses all complete cases,
%            including observations with large robust Mahalanobis distances;
%            no preliminary trimming of the radii is performed. Independent
%            directions uniform on the unit sphere are then generated, full
%            p-dimensional observations are reconstructed using the fitted
%            complete-case location and scatter, and finally the observed
%            missingness mask is imposed. Thus the observed radial tail is
%            preserved, whereas robust outlier detection and trimming are
%            reconstructed independently inside every bootstrap sample. The
%            generator does not specify a Gaussian or Student t radial family.
%            This generator is naturally paired
%            with consistencyfactor='adaptive' under an elliptical common-
%            target interpretation; it is not a general nonparametric
%            bootstrap for arbitrary skew or nonelliptical distributions.
%            Example - 'bootstraptype','empiricalradial'
%            Data Types - char | string
%
%   nsimul : Number of successful bootstrap simulations. Scalar integer.
%            The default value is 499. This option is used only when
%            bootstrap=true. If an individual bootstrap draw fails only
%            because mdTEM reports too few adaptive reference points, that
%            draw is discarded and regenerated. Other errors are rethrown.
%            Draws yielding nonfinite scalar statistics are also regenerated.
%            When omnibus calibration is active, a draw whose replicate-
%            specific discrepancy covariance or rank-selected Hausman
%            statistic is unavailable is likewise regenerated and counted in
%            out.nInvalidBootstrapOmnibus.
%            Example - 'nsimul',999
%            Data Types - double
%
%  conflev : Confidence level used to compute bootstrap confidence intervals.
%            Scalar in the interval (0,1). The default value is 0.95. This
%            option is used only when bootstrap=true.
%            Example - 'conflev',0.99
%            Data Types - double
%
%    eps0 : Positive numerical stabilizer used in the log-distance ratio.
%           The statistic is
%
%             log((d2_all + eps0) ./ (d2_cc + eps0)).
%
%           The default value is 1e-12.
%           Example - 'eps0',1e-10
%           Data Types - double
%
% aggregation : Mean aggregation rule. Character or string. Possible values are
%           'auto', 'ordinary', 'inlier' or 'fixed'. The default is 'auto'.
%
%           'auto' selects 'ordinary' when alpha=0 and 'inlier' when
%           alpha>0. Thus the default call remains the classical MDP-U
%           procedure, whereas any explicitly requested positive robustness
%           level automatically uses the recommended robust MDP-I
%           aggregation.
%
%           'ordinary' computes the two mean statistics using all complete
%           cases (MDP-U). This is the classical/theoretical aggregation
%           rule and can also be requested explicitly when alpha>0.
%
%           'inlier' computes the two mean statistics using only complete
%           cases not declared as outliers by the selected robust
%           complete-case estimator (MDP-I). This is the recommended robust
%           aggregation rule and requires alpha>0. Its primary analytical
%           calibration uses generic retained-moment directions bI,aI for
%           TDmean and their weighted analogues bLI,aLI for TLmean, propagated
%           through the joint location--scatter estimator-discrepancy sandwich.
%           The excluded fraction may converge to a positive constant under
%           the stable-inlier conditions. The scalar radial multipliers rI
%           and rLI are retained in out.asympt.inlierScaling only as the
%           retained-moment proportionality special case.
%
%           'fixed' computes the two mean statistics using only complete
%           cases whose complete-case squared Mahalanobis distance is not
%           larger than chi2inv(filterlev,p) (MDP-F). This rule is intended
%           primarily as a robust sensitivity diagnostic.
%
%           The median statistics are not affected by this option. With
%           coupledtrim=false (default), the aggregation rule acts only on
%           the final two means and does not alter the EM/TEM fit. The
%           selected mean aggregation rule is reconstructed independently
%           in every bootstrap sample.
%           Example - 'aggregation','inlier'
%           Data Types - char | string
%
% filterlev : Probability level defining the fixed robust-distance cutoff for
%           aggregation='fixed'. Scalar in the interval (0,1). The cutoff
%           is c=chi2inv(filterlev,p). The default value is 0.99.
%           This option is ignored for the other aggregation rules.
%           Example - 'filterlev',0.975
%           Data Types - double
%
% coupledtrim : Couple complete-case outlier decisions to TEM. Logical.
%           The default is false. This option is available only when
%           alpha>0 and the resolved aggregation rule is 'inlier' (for
%           example, aggregation='auto' with alpha>0, or an explicit
%           aggregation='inlier'). When true, complete cases
%           declared as outliers by the selected robust complete-case
%           estimator are forced to have zero weight in every TEM
%           concentration step. TEM still applies its own trimming rule to
%           the full data; the coupled weight is the product of the TEM
%           trimming indicator and the fixed complete-case eligibility
%           indicator. This mode is intended for sensitivity analysis and
%           is accompanied by a frozen-eligibility analytical approximation:
%           the realized eligibility mask is treated as fixed in the TEM
%           sandwich, while the variability induced by constructing that mask
%           is not differentiated. A warning is issued when this mode is used.
%           Example - 'coupledtrim',true
%           Data Types - logical
%
%  robust  : Robust estimator for the complete cases. Character, string
%            or structure. This option is used only when alpha is strictly
%            positive. If robust is a character vector or string scalar,
%            possible values are 'FS', 'MCD' or 'MM'. The default is 'FS'.
%            If robust is a structure, the following fields can be supplied:
%              robust.class = robust estimator. Possible values are 'FS',
%                             'MCD' or 'MM'. This field is required.
%              robust.eff   = nominal efficiency of the MM estimator. This
%                             field is used only if robust.class='MM'. The
%                             default value is 0.95.
%              robust.bonflev = Forward Search signal/stopping rule. This
%                             field is used only if robust.class='FS'.
%                             Use [] (default) for the standard FSM
%                             signal-detection, validation and envelope-
%                             resuperimposition rules, or 0.99 for the
%                             direct 99 percent Bonferroni-bound rule.
%            The robustness level is controlled by alpha. In the FS case,
%            floor(nComplete*(1-alpha)) is used as the FSM initial
%            monitoring step; in the MCD case alpha is passed as bdp; in
%            the MM case alpha is passed as Sbdp. For all three estimators,
%            the confidence level used to declare outliers is fixed at the
%            Bonferronized simultaneous level 0.99. In the FS case the user
%            can choose between the direct Bonferroni rule and the standard
%            envelope-resuperimposition rule through robust.bonflev.
%            Example - 'robust','MCD' | r=struct; r.class='FS'; r.bonflev=[]; 'robust',r
%            Data Types - char | string | struct
%
%      tol : Convergence tolerance passed to mdTEM. Scalar.
%            The default value is 1e-10.
%            Example - 'tol',1e-8
%            Data Types - double
%
%    plots : Flag to produce the output plot.
%            The default value is false.
%            Example - 'plots',true
%            Data Types - logical
%
% dispresults : Display commented test results. Logical scalar.
%            The default value is false. If true, mdMDPtest displays
%            the MDP procedure actually used, the principal fitting
%            options, the number of complete cases, the number entering
%            the mean statistics, and the table out.results. For classical
%            ordinary MDP-U, the separate Gaussian-information benchmark is
%            also displayed when available. If an analytical calibration
%            fails, its reason and relevant numerical-conditioning diagnostics
%            are printed rather than leaving unexplained NaNs. Small p-values
%            indicate evidence against MCAR. No automatic reject/not-reject
%            decision is printed.
%            Example - 'dispresults',true
%            Data Types - logical
%
%
%  Output:
%
%    out : Structure containing the following fields:
%
%          out.pvalueAsym  = 1 x 4 vector containing the primary asymptotic
%                            p-values in the same order as out.Tobs. The
%                            analytical theory currently concerns the two mean
%                            statistics, hence entries 1 and 3 (the medians)
%                            are NaN. For classical ordinary MDP-U, entry 2
%                            is the elliptical radial calibration of TLmean and
%                            entry 4 is the feasible general-F calibration of
%                            TDmean; Gaussian-information p-values are reported
%                            separately in out.resultsGaussian.
%          out.Tobs        = 1 x 4 vector containing the observed values of
%                            the four statistics.
%          out.results     = 4 x 6 table with row names TLmedian, TLmean,
%                            TDmedian and TDmean. The columns are Tobs,
%                            SEasym, zAsym, pvalueAsym, pvalueBoot and
%                            Calibration. The two median rows have NaN
%                            asymptotic standard errors, z values and
%                            p-values because analytical inference is
%                            currently available only for the mean summaries.
%                            The pvalueBoot column is present also when
%                            bootstrap=false, in which case it contains NaNs.
%          out.resultsGaussian = 2 x 4 table containing Tobs, SE, z and
%                            pvalue for TLmean and TDmean under the Gaussian-
%                            information benchmark. It is returned for the
%                            classical ordinary MDP-U configuration
%                            (alpha=0, aggregation='ordinary') and is empty
%                            otherwise.
%          out.omnibus     = Structure containing the Hausman-type omnibus
%                            scatter-discrepancy test. Let
%                              delta=vech(out.cov-out.covCC)
%                            and let OmegaDelta estimate the asymptotic
%                            covariance of sqrt(n)*delta. The statistic is
%                              H=n*delta'*OmegaDelta^+*delta,
%                            where ^+ denotes the Moore-Penrose inverse. Under
%                            the corresponding analytical MCAR reference,
%                            H is asymptotically chi-square with degrees of
%                            freedom equal to the statistically selected
%                            spectral rank. The rank rule retains eigenvalue j
%                            when lambda_j/lambda_1 exceeds
%                            eta_n=nEff^(-omnibusrankexp). Fields include
%                            available, stat, df, pvalue, delta, OmegaDelta,
%                            OmegaDeltaPlus, eigenvalues, relativeEigenvalues,
%                            largestEigenvalue, numericalTolerance,
%                            rankThreshold, relativeRankThreshold, rankExponent,
%                            effectiveSampleSize, effectiveSampleSizeSource,
%                            rhoLastKept, rhoFirstDiscarded, spectralGapRatio,
%                            lowerThresholdMargin, upperThresholdMargin,
%                            rankSensitivity, nullSpaceResidual, calibration
%                            and reason. The numerical tolerance is used only
%                            to diagnose roundoff-level negative eigenvalues;
%                            it does not determine the statistical rank. When
%                            bootstrap=true and the observed omnibus statistic
%                            is available, the same full fitting and covariance
%                            pipeline is rerun in every bootstrap replicate and
%                            the same spectral rank rule is recomputed. Further
%                            fields then include bootstrapAvailable,
%                            bootstrapReason, pvalueBoot, statBoot, rankBoot,
%                            relativeEigenvaluesBoot, largestEigenvalueBoot,
%                            rankThresholdBoot, relativeRankThresholdBoot,
%                            rhoLastKeptBoot, rhoFirstDiscardedBoot and
%                            rankBootFrequency.
%          out.resultsOmnibus = 1 x 5 table containing Stat, df, pvalue,
%                            pvalueBoot and Calibration for the Hausman-type
%                            omnibus test. The column pvalue is the analytical
%                            chi-square p-value; pvalueBoot is the full-pipeline
%                            bootstrap p-value. NaNs are returned when the
%                            corresponding calibration is unavailable.
%          out.bootstrap   = Value of input option bootstrap.
%          out.bootstraptype = Bootstrap null generator actually requested:
%                            'gaussian' or 'empiricalradial'.
%          out.dispresults = Value of input option dispresults.
%          out.pvalueBoot  = 1 x 4 vector containing bootstrap p-values for
%                            the four statistics. This field is present only
%                            when bootstrap=true.
%          out.Tboot       = nsimul x 4 matrix containing the bootstrap values
%                            of the four statistics. This field is present only
%                            when bootstrap=true.
%          out.alpha       = Value of input option alpha.
%          out.method      = Value of input option method.
%          out.consistencyfactor = TEM consistency-factor option actually
%                            requested. It is relevant only when alpha>0.
%          out.adaptivepool = Logical adaptive pooling option.
%          out.adaptiveminref = Minimum adaptive reference count.
%          out.omnibusrankexp = Exponent used by the omnibus spectral rank
%                            rule. The default is 1/3.
%          out.omnibusneff = Effective-sample-size specification used by the
%                            omnibus spectral rank rule.
%          out.TEMkfactor  = Scalar summary of the final TEM consistency
%                            correction. NaN for alpha=0.
%          out.TEMkinfo    = Final pattern-wise TEM consistency-factor
%                            diagnostics. Empty for alpha=0 or whenever the
%                            selected mdTEM correction does not return them.
%          out.aggregation = Mean aggregation rule actually used:
%                            'ordinary', 'inlier' or 'fixed'.
%          out.filterlev   = Value of input option filterlev.
%          out.filterCutoff = Fixed chi-square cutoff used by MDP-F. It is
%                            NaN unless aggregation='fixed'.
%          out.coupledtrim = Logical flag indicating whether estimator-
%                            declared complete-case outliers are forced to
%                            zero weight inside TEM.
%          out.coupledRows = Original row numbers forced to zero TEM weight.
%          out.nCoupledCC = Number of observed complete cases forced to zero
%                            TEM weight.
%          out.nCoupledBoot = Number of complete cases forced to zero TEM
%                            weight in each retained bootstrap sample. This
%                            field is present only when bootstrap=true.
%          out.TEMweights = Final TEM weights used for the observed data.
%          out.TEMinternalWeights = Final unconstrained TEM concentration
%                            indicators before multiplication by the fixed
%                            complete-case eligibility mask. For uncoupled
%                            fits this equals out.TEMweights.
%          out.meanWeightsCC = Logical vector of length out.nComplete. True
%                            entries identify complete cases entering the
%                            two mean statistics.
%          out.meanRows    = Original row numbers entering the two mean
%                            statistics.
%          out.nMeanCC     = Number of complete cases entering the two mean
%                            statistics.
%          out.nMeanBoot   = Number of complete cases entering the two mean
%                            statistics in each successful bootstrap sample.
%                            This field is present only when bootstrap=true.
%          out.nBootstrapAttempts = Total number of bootstrap draws attempted
%                            to obtain the requested nsimul successful draws.
%                            This field is present only when bootstrap=true.
%          out.nFailedBootstrapFits = Number of attempted bootstrap draws
%                            discarded because mdTEM returned the specific
%                            error FSDA:mdTEM:FewAdaptiveReferencePoints.
%                            Other fitting errors are not suppressed. This
%                            field is present only when bootstrap=true.
%          out.nInvalidBootstrapStats = Number of attempted bootstrap draws
%                            discarded because at least one of the four test
%                            statistics was nonfinite. This field is present
%                            only when bootstrap=true.
%          out.nInvalidBootstrapOmnibus = Number of attempted bootstrap draws
%                            regenerated because the replicate-specific
%                            scatter-discrepancy covariance or rank-selected
%                            omnibus statistic was unavailable. This field is
%                            present only when bootstrap=true.
%          out.nComplete   = Number of complete rows.
%          out.completeIdx = Logical index of complete rows.
%          out.locCC       = Complete-case estimate of location.
%          out.covCC       = Complete-case estimate of scatter.
%          out.robust      = Value of input option robust.
%          out.robustClass = Robust estimator actually used for complete
%                            cases: 'FS', 'MCD' or 'MM'.
%          out.robustEff   = MM nominal efficiency. This field is relevant
%                            only when out.robustClass='MM'.
%          out.robustBonflev = Value actually used for the FS
%                            signal/stopping rule. It is 0.99 for the
%                            direct Bonferroni rule and [] for the standard
%                            FSM envelope-resuperimposition rule. If fewer
%                            than four monitoring steps are available, []
%                            is automatically changed to 0.99.
%          out.FSrule      = Forward Search signal/stopping rule actually
%                            used: 'bonferroni' or 'resuperimposition'. It
%                            is empty if out.robustClass is not 'FS'.
%          out.outlierConflev = Simultaneous Bonferronized confidence level
%                            used to declare complete-case outliers (0.99).
%          out.outliersCC  = Complete-case units declared as outliers by
%                            the selected robust estimator.
%          out.nOutliersCC = Number of units in out.outliersCC.
%          out.d2_cc       = Mahalanobis distances computed from complete
%                            rows only.
%          out.d2_all      = Mahalanobis distances for the same complete
%                            rows when parameters are estimated from all
%                            the data through EM/TEM.
%          out.ciBoot      = 2 x 4 matrix containing the bootstrap confidence
%                            intervals for the four statistics. This field is
%                            present only when bootstrap=true.
%          out.loc          = Estimated location from EM/TEM fit on all data.
%          out.cov          = Estimated scatter from EM/TEM fit on all data.
%          out.eps0        = Numerical stabilizer used in the log-distance
%                            ratio.
%          out.asympt.inlierScaling = Stable-inlier retained-moment diagnostics
%                            when aggregation='inlier'. The primary generic
%                            directions bI,aI and bLI,aLI, together with the
%                            retained moments mI,MI,mLI,MLI, are returned.
%                            The fields rI, rLI and TLtoTDRatio remain as the
%                            proportional-moment radial special case. Additional
%                            fields quantify empirical departures from zero
%                            retained first moments and proportional retained
%                            second moments.
%          out.asympt      = Structure containing the analytical
%                            benchmark for the two selected mean statistics.
%                            For alpha=0 and aggregation='ordinary', TDmean
%                            uses the feasible general-F sandwich calibration
%                            under MCAR and finite fourth moments. Gaussianity
%                            and ellipticity are not required for this TDmean
%                            calibration. The former Gaussian-information law
%                            is retained in out.asympt.gaussian as a diagnostic
%                            specialization. Under ellipticity, TLmean is
%                            calibrated by multiplying the same general-F
%                            estimator-discrepancy variance used for TDmean
%                            by an empirical radial coefficient estimated from
%                            the complete-case Mahalanobis distances. Thus the
%                            primary TLmean calibration is elliptical but does
%                            not assume Gaussianity. The distinct Gaussian-
%                            information results are retained in
%                            out.asympt.gaussian and out.resultsGaussian.
%                            For alpha=0 and aggregation='fixed', the current
%                            analytical benchmark remains Gaussian because a
%                            general-F fixed-cutoff law has not been derived;
%                            the untrimmed general-F result is nevertheless
%                            returned in out.asympt.generalF for diagnostics.
%                            For MDP-U the ordinary first-order coefficients
%                            are used. For MDP-I, a stable estimator-induced
%                            inlier rule may exclude a nonvanishing fraction
%                            of complete cases. The primary variance uses the
%                            generic retained-moment directions bI,aI (and
%                            bLI,aLI for TLmean) with the joint location--scatter
%                            estimator-discrepancy covariance. The former rI^2
%                            and rLI^2 scalar rescalings are still returned as
%                            the proportional-moment special case, but are not
%                            used for the primary MDP-I p-values. For MDP-F the base
%                            variance is multiplied by kc^2 for TDmean and by
%                            lambda^2 for TLmean. For the Gaussian pattern
%                            branch, kc=F_{p+2}(c)/F_p(c). For the adaptive
%                            elliptical branch, kc and lambda are estimated
%                            from the robust complete-case radial distribution.
%                            Fields aggregation, filterCutoff, gammaFilter,
%                            kc, lambda, baseSigmaD2 and baseKappa document the
%                            applied scaling.
%                            out.asympt.inlierScaling is a structure with
%                            fields available, genericAvailable, nComplete,
%                            nSelected, fractionSelected, fractionExcluded,
%                            mI, MI, bI, aI, mLI, MLI, bLI, aLI, together with
%                            rI, rLI, TLtoTDRatio and radial/proportionality
%                            diagnostics. It is relevant when
%                            aggregation='inlier'. For
%                            alpha>0 and method equal to 'pri',
%                            'expScale', 'zMap', 'chiMap' or 'betaMap', the
%                            Gaussian pattern-wise TEM influence function is
%                            evaluated analytically when
%                            consistencyfactor='pattern'. For the same five
%                            methods with consistencyfactor='adaptive', the
%                            elliptical adaptive sandwich is evaluated
%                            with the empirical radial-factor influence, a
%                            kernel estimate of the radial density at each
%                            trimming boundary, and the robust complete-case
%                            scatter influence. Both pattern-specific factors
%                            (adaptivepool=false) and dimension-pooled factors
%                            (adaptivepool=true, default) are covered. In the
%                            pooled case, the influence function is formed at
%                            the complete-observation level, preserving the
%                            dependence among projections of the same row.
%                            For FS, the
%                            complete-case scatter influence is evaluated
%                            analytically conditional on the final full-sample
%                            FS classification. For MCD and MM it is currently
%                            evaluated by a delete-one jackknife. The resulting
%                            contribution is combined with the TEM influence in
%                            the end-to-end sandwich variance. With adaptive
%                            pooling, out.asympt.TEM.factorDiagnostics reports
%                            one row per pooled observed dimension, while
%                            patternDiagnostics reports the corresponding
%                            poolSize and factorGroup for each pattern.
%                            The fields TDmean and TLmean contain the
%                            asymptotic variance of sqrt(n)T (sigma2), the
%                            standard error of T (se), the studentized value
%                            (z), and the two-sided asymptotic p-value. For
%                            alpha>0, out.asympt.TEM contains diagnostics for
%                            the analytical TEM contribution and
%                            out.asympt.completeCase documents the influence
%                            evaluation, including influenceMethod and, for FS,
%                            nReferenceInliers and nReferenceOutliers. The
%                            adaptive FS branch conditions on the selected set;
%                            its data-dependent stopping/classification rule is
%                            outside the fixed-fraction FS asymptotic theorem.
%                            For alpha>0 and consistencyfactor='adaptive',
%                            out.asympt.secondOrderBenchmark contains the
%                            complete-Gaussian adaptive-TEM O(n^{-1})
%                            centering constants cMu, cQuad, cLin and cTotal.
%                            This is a theoretical benchmark only: it is not
%                            a missing-pattern second-order correction and is
%                            not used to alter analytical or bootstrap p-values.
%
%
%  More About:
%
%  Let d2_cc denote the squared Mahalanobis distances computed on the
%  complete rows using the complete-case estimates, and let d2_all denote
%  the distances for the same rows when location and scatter are estimated
%  from all the data using EM/TEM. The function monitors the following four
%  statistics:
%
%    1) median( log((d2_all + eps0) ./ (d2_cc + eps0)) );
%    2) selected mean of log((d2_all + eps0) ./ (d2_cc + eps0));
%    3) median( d2_all - d2_cc );
%    4) selected mean of d2_all - d2_cc.
%
%  With aggregation='auto', the selected mean aggregation resolves to
%  'ordinary' for alpha=0 and to 'inlier' for alpha>0. The selected means
%  use all complete cases for aggregation='ordinary', the estimator-defined
%  complete-case inliers for aggregation='inlier', and the fixed robust-
%  distance cutoff for aggregation='fixed'. The same resolved rule is
%  re-estimated in every bootstrap sample. With coupledtrim=true,
%  estimator-declared complete-case outliers are additionally forced to zero
%  weight in every TEM concentration step.
%
%  When consistencyfactor='adaptive', the complete-case observations are
%  also the reference sample for the empirical radial correction. The
%  reference location and scatter are the same complete-case estimates used
%  to compute d2_cc. Consequently the adaptive TEM correction targets the
%  complete-case scatter functional under an elliptical MCAR model. The
%  optional mask-preserving full-pipeline bootstrap can use either the fitted
%  Gaussian null generator or the null-enforcing empirical-radial elliptical
%  generator selected by bootstraptype. The latter resamples the radii of
%  all complete cases, including large robust distances, and generates fresh
%  independent spherical directions. Hence it preserves radial extremeness
%  while enforcing the elliptical angular law. In both cases the robust
%  complete-case fit, adaptive radial correction and final aggregation rule
%  are reconstructed in every bootstrap sample.
%
%  The scalar MDP statistics are directional projections of the estimator
%  discrepancy. A complementary Hausman-type omnibus diagnostic is therefore
%  also reported whenever the full scatter-discrepancy influence covariance
%  is available. With delta=vech(Sigma_all-Sigma_CC), it uses
%
%      H = n*delta'*OmegaDelta^+*delta.
%
%  Under MCAR, H has an asymptotic chi-square reference with degrees of
%  freedom equal to the statistically selected spectral rank. The relative
%  threshold is eta_n=nEff^(-omnibusrankexp), and the thresholded inverse keeps
%  only eigenvalues satisfying lambda_j/lambda_1>eta_n. Unlike the
%  one-dimensional MDP projection, this quadratic statistic is sensitive to
%  all retained first-order scatter-discrepancy directions, including
%  trace-free shape/correlation changes. With bootstrap=true, the complete
%  fitting pipeline, the discrepancy sandwich and this same rank rule are
%  recomputed independently in every retained bootstrap replicate; the
%  omnibus bootstrap p-value uses the resulting nonnegative H* statistics.
%
%  Under adaptive positive trimming, OmegaDelta is the end-to-end
%  scatter-discrepancy covariance including the adaptive radial-factor
%  influence and the robust complete-case scatter influence. Under the
%  elliptical common-target model, the same omnibus construction applies
%  to MDP-I and MDP-F: the final inlier or fixed-cutoff aggregation changes
%  the scalar MDP projection but not the underlying fitted scatter-estimator
%  discrepancy. This statement refers to the default coupledtrim=false;
%  with coupledtrim=true the eligibility rule also modifies the TEM fit and
%  the reported omnibus calibration uses the corresponding frozen-eligibility
%  working sandwich.
%
%  Small asymptotic or bootstrap p-values indicate that the change in distances
%  is larger than expected under MCAR for the corresponding calibration.
%
%  For the adaptive correction and alpha>0, the analytical sandwich assumes
%  an elliptical MCAR model and a common scatter target S0=tau*Sigma. The
%  empirical radial-factor influence and the robust complete-case scatter
%  influence are included explicitly. For the parameter-free mappings
%  'pri', 'expScale', 'zMap', 'chiMap' and 'betaMap', the analytical branch
%  covers both pattern-specific factors and dimension pooling. The mapping
%  enters through the inverse raw-distance cutoff a_g(c). In the pooled branch
%  the sampling unit remains the complete observation: projected radii from
%  the same row are not treated as independent reference observations. The
%  scatter-dependent 'detMap' adjustment is intentionally not covered by the
%  current analytical sandwich.
%
%  For alpha=0 and aggregation='ordinary', the primary analytical
%  calibration of TDmean is the feasible general-F sandwich. It estimates the
%  complete-case/all-available scatter-discrepancy influence covariance from
%  the observed data and is valid under MCAR with finite fourth moments. The
%  Gaussian-information specialization is retained in out.asympt.gaussian.
%  Under ellipticity, TLmean is a scalar multiple of TDmean to first order.
%  Its primary alpha=0 calibration therefore uses the same general-F
%  estimator-discrepancy variance multiplied by an empirical radial
%  coefficient estimated from the complete-case Mahalanobis distances.
%  This elliptical TLmean calibration is distinct from the separate Gaussian-
%  information benchmark. For alpha=0 and aggregation='fixed', the existing
%  Gaussian Tallis benchmark is retained because a general-F fixed-cutoff law
%  has not yet been derived.
%
%  For alpha>0 the analytical benchmark follows the selected mean aggregation.
%  The Gaussian pattern branch uses Tallis coefficients, whereas the adaptive
%  branch uses empirical elliptical radial coefficients. Under an elliptical
%  stable estimator-induced inlier rule, MDP-I permits a nonvanishing
%  excluded fraction. The generic calibration uses retained first and second
%  moments directly, so a non-centrally-symmetric retained distribution may
%  contribute a first-order location term. When retained first moments vanish
%  and retained second moments are proportional to the common target, this
%  reduces to the scalar rI^2 and rLI^2 radial rescalings. MDP-F
%  multiplies the base
%  mean-difference variance by kc^2 and the base stabilized-log-ratio variance
%  by the appropriate truncated radial coefficient lambda^2. When
%  coupledtrim=true, mdMDPtest still reports an
%  analytical p-value based on a frozen-eligibility working sandwich: the
%  realized complete-case eligibility mask is treated as fixed, and its
%  data-dependent construction is not differentiated. A warning documents
%  this approximation. Bootstrap calibration can be added with bootstrap=true.
%
% See also: mdEM, mdTEM, mdImputeCondMean.m, mdPartialMD.m, mdPartialMD2full
%
% References:
%
% Little, R. J. A., & Rubin, D. B. (2019). Statistical Analysis with
% Missing Data (3rd ed.). Hoboken, NJ: John Wiley & Sons.
% van Buuren, S. (2018). Flexible Imputation of Missing Data (2nd ed.).
% Boca Raton, FL: Chapman & Hall/CRC (Taylor & Francis Group).
% Templ, M. (2023). Visualization and Imputation of Missing Values: With
% Applications in R. Cham, Switzerland: Springer Nature.
% Hausman, J. A. (1978). Specification Tests in Econometrics. Econometrica,
% 46, 1251-1271.
%
%
% Copyright 2008-2026.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mdMDPtest')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example 1: Classical MDP-U with default options.
    % The numerical default alpha=0 is unchanged. With the default
    % aggregation='auto', alpha=0 resolves to aggregation='ordinary'.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X);
    % Display the two statistics TLmean and TDmean 
    % together with their asymptotic p-values.
    disp(out.results(["TLmean" "TDmean"],{'Tobs','SEasym','zAsym','pvalueAsym','Calibration'}))
%}

%{
    %% Example 2: Recommended robust MDP-I analysis.
    % A positive alpha activates robust complete-case/TEM estimation. With
    % aggregation='auto' (default), the two mean statistics automatically
    % use the estimator-defined complete-case inliers (MDP-I). The default
    % TEM correction is distribution-adaptive with pooling by observed
    % dimension. The value alpha=0.25 is illustrative and is not a new
    % numerical default.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25);
    % Display the observed statistics and asymptotic p-values with row names.
    disp(out.results(:,{'Tobs','SEasym','zAsym','pvalueAsym','Calibration'}))
    % Display the generic retained-moment directions and radial special-case diagnostics.
    disp(out.asympt.inlierScaling)
%}

%{
    %% Example 3: Simulated data under MCAR.
    % Generate Gaussian data with MCAR missingness and apply the test.
    rng(1)
    n = 300;
    p = 5;
    rho = 0.5;
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    mu = zeros(1,p);

    Yfull = mvnrnd(mu,Sigma,n);

    missRate = 0.10;
    missMask = rand(n,p) < missRate;
    Y = Yfull;
    Y(missMask) = NaN;
    % Show also the output plot and request bootstrap calibration.
    out = mdMDPtest(Y,'bootstrap',true,'nsimul',199,'plots',true);

    % The table identifies each statistic and reports both calibrations.
    disp(out.results)
%}


%{
    %% Example 4: Comparison of several trimming levels.
    % Compare bootstrap p-values for the four named MDP statistics.
    rng(1)
    n = 300;
    p = 5;
    rho = 0.5;
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    mu = zeros(1,p);
    Yfull = mvnrnd(mu,Sigma,n);
    missRate = 0.10;
    missMask = rand(n,p) < missRate;
    Y = Yfull;
    Y(missMask) = NaN;    
    Alpha = [0 0.10 0.25 0.50]';
    pval = zeros(length(Alpha),4);
    for i=1:length(Alpha)
        out = mdMDPtest(Y,'alpha',Alpha(i),'bootstrap',true,'nsimul',199);
        pval(i,:) = out.results.pvalueBoot';
    end
    pvalTable = table(Alpha,pval(:,1),pval(:,2),pval(:,3),pval(:,4), ...
        'VariableNames',{'alpha','pBoot_TLmedian','pBoot_TLmean', ...
        'pBoot_TDmedian','pBoot_TDmean'});
    disp(pvalTable)
%}

%{
    %% Example 5: Different rescaling method.
    % Use method betaMap instead of the default pri.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'method','betaMap', ...
        'bootstrap',true,'nsimul',199);
    % Display statistics, asymptotic p-values and bootstrap p-values together.
    disp(out.results)
%}

%{
    % Example 6: MM estimator for the complete cases.
    % Use an MM estimator with 95 percent nominal efficiency. The
    % preliminary S estimator has breakdown point alpha=0.25. Outliers are
    % declared using the fixed Bonferronized simultaneous level 0.99.
    load cows2026
    X = cows2026{:,:};
    r = struct;
    r.class = 'MM';
    r.eff = 0.95;
    out = mdMDPtest(X,'alpha',0.25,'robust',r);
    % Display the four named statistics and their asymptotic p-values.
    disp(out.results(:,{'Tobs','SEasym','zAsym','pvalueAsym','Calibration'}))
%}

%{
    % Example 7: Standard FSM envelope-resuperimposition rule.
    % Use the standard FSM signal-detection and validation procedure rather
    % than the direct Bonferroni-bound stopping rule.
    load cows2026
    X = cows2026{:,:};
    r = struct;
    r.class = 'FS';
    r.bonflev = [];
    out = mdMDPtest(X,'alpha',0.25,'robust',r);
    fprintf('Forward Search rule used: %s\n',out.FSrule)
    disp(out.results(:,{'Tobs','SEasym','zAsym','pvalueAsym','Calibration'}))
%}

%{
    % Example 8: Explicit estimator-defined inlier mean (MDP-I).
    % The medians use all complete cases. The two means use only complete
    % cases not declared as outliers by the MCD complete-case fit.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'robust','MCD', ...
        'aggregation','inlier');
    fprintf('Complete cases entering the mean statistics: %d\n',out.nMeanCC)
    disp(out.results(:,{'Tobs','SEasym','zAsym','pvalueAsym','Calibration'}))
    % rI and rLI allow the excluded fraction to remain nonvanishing.
    disp(out.asympt.inlierScaling)
%}

%{
    % Example 9: Fixed robust-distance sensitivity diagnostic (MDP-F).
    % Use the chi-square 0.99 cutoff on the complete-case robust distances.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'robust','MCD', ...
        'aggregation','fixed','filterlev',0.99);
    fprintf('Fixed robust-distance cutoff: %g\n',out.filterCutoff)
    fprintf('Complete cases entering the mean statistics: %d\n',out.nMeanCC)
    disp(out.results(:,{'Tobs','SEasym','zAsym','pvalueAsym','Calibration'}))

    % Display the MDP-F first-order coefficients with explicit names.
    coefTable = table(out.asympt.kc,out.asympt.lambda, ...
        'VariableNames',{'kc_TDmean','lambda_TLmean'});
    disp(coefTable)

    % Display the two asymptotic variances with explicit test names.
    varTable = table(out.asympt.TDmean.sigma2,out.asympt.TLmean.sigma2, ...
        'VariableNames',{'sigma2_TDmean','sigma2_TLmean'});
    disp(varTable)
%}

%{
    % Example 10: Coupled complete-case outlier trimming in TEM.
    % This sensitivity mode uses the MDP-I aggregation rule and additionally
    % forces estimator-declared complete-case outliers to zero TEM weight.
    % Coupled trimming currently requires the Gaussian pattern correction,
    % which is therefore requested explicitly below.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'robust','MCD', ...
        'aggregation','inlier','consistencyfactor','pattern', ...
        'coupledtrim',true,'bootstrap',true,'nsimul',199);
    fprintf('Complete cases forced to zero TEM weight: %d\n',out.nCoupledCC)
    disp('Rows forced to zero TEM weight:')
    disp(out.coupledRows)
    % Display statistics and both p-value calibrations with explicit row names.
    disp(out.results)
%}

%{
    % Example 11: Distribution-adaptive pooled TEM correction.
    % Use the MCD complete-case fit as the reference geometry for empirical
    % radial consistency factors. With adaptivepool=true (default), projected
    % reference radii are pooled across patterns of the same observed
    % dimension. The same pipeline is rerun independently inside every
    % bootstrap sample.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'robust','MCD', ...
        'consistencyfactor','adaptive');
    disp(out.TEMkinfo)
    disp(out.results(:,{'Tobs','SEasym','zAsym','pvalueAsym','Calibration'}))
    disp(out.asympt.TDmean)
    disp(out.asympt.TEM.patternDiagnostics)
%}


%{
    %% Example 12: General-F analytical calibration for untrimmed MDP-U.
    % For alpha=0 and ordinary aggregation, TDmean uses the feasible
    % general-F sandwich as the primary analytical calibration. The
    % Gaussian-information specialization is retained as a diagnostic.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0,'aggregation','ordinary');
    % Display the primary analytical results. TLmean uses the elliptical
    % radial extension of the general-F sandwich; TDmean uses general-F.
    disp(out.results)
    % Display separately the Gaussian-information benchmark.
    disp(out.resultsGaussian)
%}

%{
    %% Example 13: Display a commented summary and the results table.
    % With dispresults=true, the procedure, fitting options, sample sizes,
    % and the named results table are displayed automatically.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'dispresults',true,'bootstrap',true);
    % out.results contains the primary calibration and out.resultsGaussian
    % contains the distinct Gaussian benchmark when alpha=0 and MDP-U is used.
%}

%{
    %% Example 14: Empirical-radial elliptical bootstrap.
    % Calibrate robust MDP-I by resampling the empirical complete-case radial
    % law rather than imposing a Gaussian radial distribution. The observed
    % missingness mask is preserved and the full robust MDP pipeline is rerun
    % independently in every bootstrap replicate.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'bootstrap',true, ...
        'bootstraptype','empiricalradial','nsimul',999);
    fprintf('Bootstrap generator: %s\n',out.bootstraptype)
    disp(out.results)
%}

%{
    %% Example 15: Hausman-type omnibus scatter-discrepancy diagnostic.
    % Compare the complete-case and all-available scatter estimators through
    % the full quadratic estimator discrepancy rather than through only the
    % one-dimensional MDP projection. With alpha>0 the default mean
    % aggregation is MDP-I, but the final inlier aggregation does not alter
    % the omnibus statistic because the latter acts directly on the fitted
    % scatter-estimator discrepancy.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'bootstrap',true,'nsimul',199);
    disp(out.resultsOmnibus)
    % Detailed observed and bootstrap rank diagnostics are stored in out.omnibus.
    disp(out.omnibus.rankSensitivity)
    if out.omnibus.bootstrapAvailable
        disp(out.omnibus.rankBootFrequency)
    end
%}

%% Beginning of code

if istable(Y)
    Y=Y{:,:};
end


if ~ismatrix(Y) || ~isnumeric(Y)
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Input argument Y must be a numeric matrix.');
end

% Default options
options = struct;
options.alpha   = 0;
options.method  = 'pri';
options.bootstrap = false;
options.bootstraptype = 'gaussian';
options.nsimul  = 499;
options.conflev = 0.95;
options.tol     = 1e-10;
options.plots   = false;
options.dispresults = false;
options.robust  ="FS";
options.eps0    = 1e-12;
options.aggregation = 'auto';
options.filterlev = 0.99;
options.coupledtrim = false;
options.consistencyfactor = 'adaptive';
options.adaptivepool = true;
options.adaptiveminref = 20;
options.omnibusrankexp = 1/3;
options.omnibusneff = 'n';

% Check supplied options
if ~isempty(varargin)
    if mod(length(varargin),2) ~= 0
        error('FSDA:mdMDPtest:WrongInputOpt', ...
            'Optional arguments must be supplied in name/value pairs.');
    end

    UserOptions = varargin(1:2:end);
    if ~isempty(UserOptions)
        aux.chkoptions(options,UserOptions);
    end

    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
end

alpha   = options.alpha;
method  = char(options.method);
bootstrap = options.bootstrap;
bootstraptype = local_parse_bootstraptype(options.bootstraptype);
nsimul  = options.nsimul;
conflev = options.conflev;
tol     = options.tol;
plots   = options.plots;
dispresults = options.dispresults;
robust  = options.robust;
eps0 = options.eps0;
aggregation = local_parse_aggregation(options.aggregation);
filterlev = options.filterlev;
coupledtrim = options.coupledtrim;
consistencyfactor = local_parse_consistencyfactor(options.consistencyfactor);
adaptivepool = options.adaptivepool;
adaptiveminref = options.adaptiveminref;
omnibusrankexp = options.omnibusrankexp;
omnibusneff = options.omnibusneff;

if ~isscalar(alpha) || ~isnumeric(alpha) || alpha < 0 || alpha > 0.5
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option alpha must be a scalar in the interval [0,0.5].');
end

% Resolve the automatic aggregation policy only after alpha has been
% validated. The numerical default alpha=0 is deliberately unchanged:
%   alpha=0  -> classical/theoretical MDP-U ('ordinary');
%   alpha>0  -> recommended robust MDP-I ('inlier').
% Explicit 'ordinary', 'inlier' and 'fixed' requests are left unchanged.
if strcmp(aggregation,'auto')
    if alpha > 0
        aggregation = 'inlier';
    else
        aggregation = 'ordinary';
    end
end

if ~(islogical(bootstrap) && isscalar(bootstrap))
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option bootstrap must be a logical scalar.');
end

if ~(islogical(dispresults) && isscalar(dispresults))
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option dispresults must be a logical scalar.');
end

% nsimul is deliberately ignored when bootstrap=false. This permits the
% cheap asymptotic-only route without imposing bootstrap-specific checks.
if bootstrap && (~isscalar(nsimul) || ~isnumeric(nsimul) || ...
        ~isfinite(nsimul) || nsimul <= 0 || nsimul ~= floor(nsimul))
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option nsimul must be a positive integer when bootstrap=true.');
end

if ~isscalar(conflev) || conflev <= 0 || conflev >= 1
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option conflev must be a scalar in the interval (0,1).');
end

if ~isscalar(tol) || tol <= 0
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option tol must be a positive scalar.');
end

if ~isscalar(eps0) || ~isnumeric(eps0) || ~isfinite(eps0) || eps0 <= 0
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option eps0 must be a finite positive scalar.');
end

if ~isscalar(filterlev) || ~isnumeric(filterlev) || ...
        ~isfinite(filterlev) || filterlev <= 0 || filterlev >= 1
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option filterlev must be a finite scalar in the interval (0,1).');
end

if ~(islogical(coupledtrim) && isscalar(coupledtrim))
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option coupledtrim must be a logical scalar.');
end

if ~(islogical(adaptivepool) && isscalar(adaptivepool))
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option adaptivepool must be a logical scalar.');
end

if ~isscalar(adaptiveminref) || ~isnumeric(adaptiveminref) || ...
        ~isfinite(adaptiveminref) || adaptiveminref < 1 || ...
        adaptiveminref ~= floor(adaptiveminref)
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option adaptiveminref must be a positive integer.');
end

if ~isscalar(omnibusrankexp) || ~isnumeric(omnibusrankexp) || ...
        ~isfinite(omnibusrankexp) || omnibusrankexp <= 0 || ...
        omnibusrankexp >= 0.5
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option omnibusrankexp must be a finite scalar in the open interval (0,0.5).');
end

if coupledtrim && ~strcmp(consistencyfactor,'pattern')
    error('FSDA:mdMDPtest:CoupledConsistencyFactor', ...
        ['The current coupledtrim sensitivity implementation supports only ' ...
        'consistencyfactor=''pattern''.']);
end

if coupledtrim && alpha == 0
    error('FSDA:mdMDPtest:CoupledNeedsRobust', ...
        'Option coupledtrim=true requires alpha>0.');
end

if coupledtrim && ~strcmp(aggregation,'inlier')
    error('FSDA:mdMDPtest:CoupledNeedsInlier', ...
        ['Option coupledtrim=true is currently defined only for ' ...
        'aggregation=''inlier''.']);
end

if strcmp(aggregation,'inlier') && alpha == 0
    error('FSDA:mdMDPtest:InlierNeedsRobust', ...
        ['Aggregation ''inlier'' requires alpha>0 because the complete-case ' ...
        'outlier set must be supplied by a robust complete-case fit.']);
end

% Parse the complete-case robust estimator. A simple character/string
% specification is retained for convenience; a structure permits
% estimator-specific tuning such as the MM efficiency.
[robustClass,robustEff,robustBonflev] = local_parse_robust(robust);

[n,p] = size(Y);
maskMiss = isnan(Y);
completeIdx = all(~maskMiss,2);
nComplete = sum(completeIdx);

% Resolve the effective sample size used only by the omnibus spectral rank
% threshold. This does not alter fitting or the discrepancy covariance.
[nEffOmnibus,nEffOmnibusSource,omnibusneff] = ...
    local_resolve_omnibus_neff(omnibusneff,n,nComplete);

if nComplete < p + 2
    error('FSDA:mdMDPtest:TooFewCompleteRows', ...
        ['Too few complete rows to compute the reference complete-case ' ...
        'covariance matrix.']);
end

if alpha > 0 && strcmpi(robustClass,'FS')

    h = floor(nComplete*(1-alpha));

    if h < p+1
        error('FSDA:mdMDPtest:TooSmallFSInit', ...
            ['The value of alpha produces an FSM initial monitoring step ' ...
            'smaller than p+1. Decrease alpha or increase the number ' ...
            'of complete observations.']);
    end

    nMon = nComplete-h;

    if nMon < 3
        error('FSDA:mdMDPtest:TooLargeFSInit', ...
            ['The value of alpha leaves fewer than three Forward Search ' ...
            'monitoring steps. Increase alpha.']);

    elseif nMon < 4 && isempty(robustBonflev)

        robustBonflev = 0.99;

        warning('FSDA:mdMDPtest:FSBonferroniFallback', ...
            ['Fewer than four Forward Search monitoring steps are ' ...
            'available. The direct 99 percent Bonferroni stopping ' ...
            'rule is used instead of the standard FSM ' ...
            'envelope-resuperimposition rule.']);
    end
end


% Extraction of complete cases
Ycc = Y(completeIdx,:);

% Robust mu and Sigma based on complete cases
[muCC,SigCC,outRob] = local_complete_case_fit(Ycc,alpha, ...
    robustClass,robustEff,robustBonflev);

% Distances based on complete rows only
d2_cc = mahalFS(Ycc, muCC, SigCC);

% Construct the complete-case eligibility mask used by the optional
% coupled TEM sensitivity mode. With coupledtrim=false no row is forced out.
forcedZeroTEM = false(n,1);
if coupledtrim
    ccOutlierMask = local_cc_outlier_mask(outRob,nComplete);
    completeRows = find(completeIdx);
    forcedZeroTEM(completeRows(ccOutlierMask)) = true;
end

% Distances based on EM/TEM fit using all rows. The fitted object is kept
% because the alpha>0 analytical benchmark needs the final TEM weights and
% pattern-wise trimming information. In coupled mode the fixed eligibility
% mask is multiplied by TEM's own concentration weights at every iteration.
[d2_all_cc, muHat, SigHat, outFit] = local_fit_and_get_complete_distances( ...
    Y, completeIdx, alpha, method, tol, forcedZeroTEM, ...
    consistencyfactor, Ycc, muCC, SigCC, adaptivepool, adaptiveminref);

% Construct the observed-data aggregation rule for the two mean statistics.
% The median statistics always use all complete cases.
[meanWeightsCC,aggInfo] = local_aggregation_weights(d2_cc,outRob,alpha, ...
    aggregation,filterlev,p);

if aggInfo.nSelected == 0
    error('FSDA:mdMDPtest:NoAggregationRows', ...
        'The selected mean aggregation rule retains no complete cases.');
end

% Observed statistics
Tobs = local_statistic(d2_cc,d2_all_cc,eps0,meanWeightsCC);

% Keep the untrimmed ordinary TDmean available even when aggregation='fixed'.
% The feasible general-F theorem applies to this MDP-U statistic only.
TDordinary = mean(d2_all_cc-d2_cc);

% Analytical benchmark for the selected mean statistics.
%
% For alpha=0 and ordinary aggregation, TDmean uses the feasible general-F
% sandwich as the primary analytical calibration. The Gaussian information
% specialization is retained as a diagnostic. For alpha=0 with a fixed cutoff,
% the present analytical law remains Gaussian because a general-F fixed-cutoff
% result has not yet been derived. For alpha>0 the existing Gaussian-pattern
% or adaptive-elliptical TEM branches are used.
if coupledtrim
    % Frozen-eligibility working sandwich. The realized externally imposed
    % complete-case eligibility mask is held fixed. Its effect on the TEM
    % estimating equation is represented through the eligible pattern mass
    % and final product weights, but the variability induced by constructing
    % the eligibility mask is deliberately not differentiated.
    asympt = local_tem_asymptotic(Y, maskMiss, completeIdx, outFit, ...
        alpha, method, robustClass, robustEff, robustBonflev, ...
        Tobs, eps0, outRob, ~forcedZeroTEM);
    asympt = local_apply_aggregation_asymptotic(asympt,aggregation,aggInfo, ...
        p,eps0,n,Tobs,d2_cc,Ycc,muCC,SigCC);
    asympt.mode = ['TEM influence-function sandwich with generic stable-inlier ' ...
        'retained-moment calibration (frozen eligibility)'];
    asympt.theoryStatus = [asympt.theoryStatus ' With coupledtrim=true, ' ...
        'the realized complete-case eligibility mask is treated as fixed; ' ...
        'the variability induced by its data-dependent construction is not ' ...
        'included in the sandwich approximation.'];
    warning('FSDA:mdMDPtest:CoupledAsymptoticApproximation', ...
        ['With coupledtrim=true the analytical p-value uses a frozen-' ...
        'eligibility approximation. The realized complete-case eligibility ' ...
        'mask is treated as fixed and its selection variability is not ' ...
        'included in the current sandwich calculation.']);
else
    if alpha == 0
        preferGeneralF = strcmp(aggregation,'ordinary');
        asympt = local_classical_asymptotic(Y, maskMiss, completeIdx, ...
            muHat, SigHat, nComplete, Tobs, TDordinary, eps0, d2_cc, ...
            preferGeneralF);

        % The feasible general-F result currently applies only to ordinary
        % untrimmed TDmean. For the fixed cutoff, retain the Gaussian Tallis
        % scaling of the former analytical benchmark.
        if ~preferGeneralF
            asympt = local_apply_aggregation_asymptotic(asympt,aggregation, ...
                aggInfo,p,eps0,n,Tobs,d2_cc,Ycc,muCC,SigCC);
        end
    elseif strcmp(consistencyfactor,'pattern')
        asympt = local_tem_asymptotic(Y, maskMiss, completeIdx, outFit, ...
            alpha, method, robustClass, robustEff, robustBonflev, ...
            Tobs, eps0, outRob, []);
        asympt = local_apply_aggregation_asymptotic(asympt,aggregation,aggInfo, ...
            p,eps0,n,Tobs,d2_cc,Ycc,muCC,SigCC);
    elseif strcmp(consistencyfactor,'adaptive')
        asympt = local_tem_adaptive_asymptotic(Y, maskMiss, completeIdx, ...
            outFit, alpha, method, robustClass, robustEff, robustBonflev, ...
            Tobs, eps0, muCC, SigCC, adaptivepool, adaptiveminref, outRob);
        asympt = local_apply_aggregation_asymptotic(asympt,aggregation,aggInfo, ...
            p,eps0,n,Tobs,d2_cc,Ycc,muCC,SigCC);
    else
        asympt = local_unavailable_asymptotic(nComplete/n);
        asympt.reason = ['For alpha>0 the analytical TEM benchmark is ' ...
            'currently implemented only for consistencyfactor=''pattern'' ' ...
            'or ''adaptive''.'];
        asympt.mode = 'TEM bootstrap calibration';
        asympt.theoryStatus = 'Bootstrap calibration for the selected TEM correction.';
        asympt = local_apply_aggregation_asymptotic(asympt,aggregation,aggInfo, ...
            p,eps0,n,Tobs,d2_cc,Ycc,muCC,SigCC);
    end
end

% Complete-Gaussian adaptive-TEM second-order benchmark.  This is exposed
% only as a diagnostic theoretical reference.  It is deliberately not used
% to recenter the observed statistic or the analytical/bootstrap p-values,
% because the closed-form cLin result has not yet been extended to genuine
% missing-pattern mixtures or to the robust complete-case reference fit.
if alpha > 0 && strcmp(consistencyfactor,'adaptive')
    asympt.secondOrderBenchmark = mdMDPsecondOrder(p,alpha);
else
    asympt.secondOrderBenchmark = [];
end

% Hausman-type omnibus estimator-discrepancy diagnostic. This uses the full
% asymptotic covariance of vech(Sigma_all-Sigma_CC), when available, rather
% than the one-dimensional MDP projection. The statistic is independent of
% the final mean aggregation rule because it compares the two fitted scatter
% estimators themselves.
omnibus = local_omnibus_test(asympt,SigHat,SigCC,n,alpha, ...
    consistencyfactor,coupledtrim,omnibusrankexp,nEffOmnibus, ...
    nEffOmnibusSource);

% Bootstrap-specific omnibus fields are initialized even when bootstrap is
% not requested so that the omnibus output has a stable structure.
omnibus.bootstrapAvailable = false;
omnibus.bootstrapReason = 'Bootstrap calibration not requested.';
omnibus.bootstrapCalibration = '';
omnibus.pvalueBoot = NaN;
omnibus.statBoot = [];
omnibus.rankBoot = [];
omnibus.relativeEigenvaluesBoot = [];
omnibus.largestEigenvalueBoot = [];
omnibus.rankThresholdBoot = [];
omnibus.relativeRankThresholdBoot = [];
omnibus.rhoLastKeptBoot = [];
omnibus.rhoFirstDiscardedBoot = [];
omnibus.rankBootFrequency = table();

% Assemble the cheap analytical p-value vector in the same order as Tobs.
% The present analytical theory concerns the two mean statistics only.
pvalueAsym = NaN(1,4);
if isfield(asympt,'TLmean') && isstruct(asympt.TLmean) && ...
        isfield(asympt.TLmean,'pvalue')
    pvalueAsym(2) = asympt.TLmean.pvalue;
end
if isfield(asympt,'TDmean') && isstruct(asympt.TDmean) && ...
        isfield(asympt.TDmean,'pvalue')
    pvalueAsym(4) = asympt.TDmean.pvalue;
end

% Optional mask-preserving full-pipeline bootstrap under MCAR. The
% expensive resampling pipeline is entered only when bootstrap=true. The
% generator is selected through bootstraptype; all subsequent fitting and
% aggregation steps are identical for the two generators. When the observed
% omnibus statistic is available, each retained replicate also reconstructs
% the full scatter-discrepancy sandwich and reapplies the same rank-selection
% rule before its bootstrap Hausman statistic is stored.
Tboot = [];
nMeanBoot = [];
nCoupledBoot = [];
pvalueBoot = [];
ciBoot = [];
Hboot = [];
rankBoot = [];
relativeEigenvaluesBoot = [];
largestEigenvalueBoot = [];
rankThresholdBoot = [];
relativeRankThresholdBoot = [];
rhoLastKeptBoot = [];
rhoFirstDiscardedBoot = [];

if bootstrap
    Tboot = NaN(nsimul,4);
    nMeanBoot = NaN(nsimul,1);
    nCoupledBoot = NaN(nsimul,1);
    if omnibus.available
        Hboot = NaN(nsimul,1);
        rankBoot = NaN(nsimul,1);
        relativeEigenvaluesBoot = NaN(nsimul,omnibus.dimension);
        largestEigenvalueBoot = NaN(nsimul,1);
        rankThresholdBoot = NaN(nsimul,1);
        relativeRankThresholdBoot = NaN(nsimul,1);
        rhoLastKeptBoot = NaN(nsimul,1);
        rhoFirstDiscardedBoot = NaN(nsimul,1);
        omnibus.bootstrapReason = '';
    else
        omnibus.bootstrapReason = ['Observed omnibus statistic unavailable: ' ...
            omnibus.reason];
    end

    % nsimul denotes the number of successful bootstrap replicates. A draw
    % that fails only because mdTEM cannot estimate an adaptive factor from
    % enough reference radii is discarded and regenerated. Other errors are
    % deliberately rethrown. Nonfinite statistic draws are also regenerated.
    % The safety margin prevents an infinite loop in pathological cases.
    nBootstrapAttempts = 0;
    nFailedBootstrapFits = 0;
    nInvalidBootstrapStats = 0;
    nInvalidBootstrapOmnibus = 0;
    maxBootstrapAttempts = nsimul + max(100,ceil(0.10*nsimul));

    % Fitted complete-case geometry used by both null generators. Rgen is
    % upper triangular and satisfies Rgen'*Rgen=SigGen. Since observations
    % are stored as row vectors, multiplying a spherical row vector by Rgen
    % gives the required fitted scatter geometry.
    SigGen = local_make_spd(SigCC);
    Rgen = chol(SigGen,'upper');

    % The empirical-radial generator keeps the observed complete-case radial
    % law but replaces the empirical angular distribution by independent
    % directions uniform on the unit sphere. All complete cases contribute
    % to the empirical radial distribution; robust selection is deliberately
    % reconstructed afresh inside each bootstrap replicate.
    QrefBoot = [];
    nQBoot = 0;
    if strcmp(bootstraptype,'empiricalradial')
        QrefBoot = mahalFS(Ycc,muCC,SigGen);
        QrefBoot = QrefBoot(:);
        if isempty(QrefBoot) || any(~isfinite(QrefBoot)) || any(QrefBoot < 0)
            error('FSDA:mdMDPtest:InvalidRadialReference', ...
                ['Unable to construct the empirical-radial bootstrap because ' ...
                'the complete-case squared distances are empty, negative or nonfinite.']);
        end
        nQBoot = numel(QrefBoot);
    end

    j = 0;
    while j < nsimul
        if nBootstrapAttempts >= maxBootstrapAttempts
            error('FSDA:mdMDPtest:TooManyBootstrapRetries', ...
                ['Unable to obtain %d valid bootstrap replicates within %d ' ...
                'attempts. Adaptive-TEM scarcity failures: %d; draws with ' ...
                'nonfinite statistics: %d; unavailable omnibus draws: %d.'], ...
                nsimul,maxBootstrapAttempts,nFailedBootstrapFits, ...
                nInvalidBootstrapStats,nInvalidBootstrapOmnibus);
        end
        nBootstrapAttempts = nBootstrapAttempts + 1;

        try
            % Generate latent full data under the selected MCAR null generator.
            switch bootstraptype
                case 'gaussian'
                    YfullStar = randn(n,p) * Rgen + muCC(:)';

                case 'empiricalradial'
                    indQ = randi(nQBoot,n,1);
                    qStar = QrefBoot(indQ);
                    Udir = local_random_spherical_directions(n,p);
                    Xsphere = Udir.*sqrt(qStar);
                    YfullStar = Xsphere * Rgen + muCC(:)';
            end

            % Impose exactly the observed missingness pattern.
            Ystar = YfullStar;
            Ystar(maskMiss) = NaN;

            % Complete-case reference distances in bootstrap world.
            YccStar = YfullStar(completeIdx,:);

            [muCCStar,SigCCStar,outRobStar] = local_complete_case_fit(YccStar,alpha, ...
                robustClass,robustEff,robustBonflev);

            d2_cc_star = mahalFS(YccStar, muCCStar, SigCCStar);

            % Reconstruct the coupled complete-case eligibility mask independently
            % in every bootstrap sample.
            forcedZeroStar = false(n,1);
            if coupledtrim
                ccOutlierMaskStar = local_cc_outlier_mask(outRobStar,nComplete);
                completeRows = find(completeIdx);
                forcedZeroStar(completeRows(ccOutlierMaskStar)) = true;
            end
            nCoupledStar = sum(forcedZeroStar);

            % EM/TEM distances and fitted geometry for the same complete rows.
            % The fitted object is retained because the omnibus bootstrap must
            % reconstruct the replicate-specific discrepancy sandwich.
            [d2_all_cc_star,muHatStar,SigHatStar,outFitStar] = ...
                local_fit_and_get_complete_distances( ...
                Ystar, completeIdx, alpha, method, tol, forcedZeroStar, ...
                consistencyfactor, YccStar, muCCStar, SigCCStar, ...
                adaptivepool, adaptiveminref);

            % Reconstruct the selected aggregation rule inside this bootstrap sample.
            [meanWeightsStar,aggInfoStar] = local_aggregation_weights( ...
                d2_cc_star,outRobStar,alpha,aggregation,filterlev,p);
            nMeanStar = aggInfoStar.nSelected;

            % Compute the four statistics. A draw with any nonfinite statistic
            % is discarded below and regenerated, so successful output always
            % contains exactly nsimul valid rows.
            Tstar = local_statistic(d2_cc_star,d2_all_cc_star,eps0, ...
                meanWeightsStar);

            % Full-pipeline omnibus calibration. When the observed omnibus
            % statistic is available, reconstruct the replicate-specific
            % scatter-discrepancy covariance and apply exactly the same
            % statistical spectral rank rule used for the observed sample.
            if omnibus.available
                asymptStarOmnibus = local_bootstrap_omnibus_asymptotic( ...
                    Ystar,maskMiss,completeIdx,outFitStar,alpha,method, ...
                    robustClass,robustEff,robustBonflev,Tstar,eps0, ...
                    outRobStar,coupledtrim,forcedZeroStar,consistencyfactor, ...
                    muCCStar,SigCCStar,adaptivepool,adaptiveminref, ...
                    d2_cc_star,d2_all_cc_star,muHatStar,SigHatStar);

                [nEffOmnibusStar,nEffOmnibusSourceStar] = ...
                    local_resolve_omnibus_neff(omnibusneff,n,sum(completeIdx));
                omnibusStar = local_omnibus_test(asymptStarOmnibus, ...
                    SigHatStar,SigCCStar,n,alpha,consistencyfactor, ...
                    coupledtrim,omnibusrankexp,nEffOmnibusStar, ...
                    nEffOmnibusSourceStar);
            else
                omnibusStar = [];
            end

        catch ME
            % This is the one expected numerical scarcity failure identified
            % in the adaptive-TEM validation study. It concerns an individual
            % bootstrap draw, not the observed-data fit. Regenerate the draw.
            if strcmp(ME.identifier,'FSDA:mdTEM:FewAdaptiveReferencePoints')
                nFailedBootstrapFits = nFailedBootstrapFits + 1;
                continue
            end
            rethrow(ME)
        end

        if any(~isfinite(Tstar))
            nInvalidBootstrapStats = nInvalidBootstrapStats + 1;
            continue
        end

        if omnibus.available && (~isstruct(omnibusStar) || ...
                ~omnibusStar.available || ~isfinite(omnibusStar.stat) || ...
                omnibusStar.stat < 0)
            nInvalidBootstrapOmnibus = nInvalidBootstrapOmnibus + 1;
            continue
        end

        % Retain only successful draws. The storage index is incremented only
        % here, so Tboot, nMeanBoot and nCoupledBoot always have nsimul rows.
        j = j + 1;
        Tboot(j,:) = Tstar;
        nMeanBoot(j) = nMeanStar;
        nCoupledBoot(j) = nCoupledStar;
        if omnibus.available
            Hboot(j) = omnibusStar.stat;
            rankBoot(j) = omnibusStar.rank;
            relativeEigenvaluesBoot(j,:) = omnibusStar.relativeEigenvalues(:)';
            largestEigenvalueBoot(j) = omnibusStar.largestEigenvalue;
            rankThresholdBoot(j) = omnibusStar.rankThreshold;
            relativeRankThresholdBoot(j) = omnibusStar.relativeRankThreshold;
            rhoLastKeptBoot(j) = omnibusStar.rhoLastKept;
            rhoFirstDiscardedBoot(j) = omnibusStar.rhoFirstDiscarded;
        end
    end

    % Bootstrap p-values for the four scalar statistics.
    pvalueBoot = (1 + sum(abs(Tboot) >= abs(Tobs),1)) / (size(Tboot,1) + 1);

    % Rank-selected omnibus bootstrap p-value. H is nonnegative, so the
    % bootstrap comparison uses its upper tail rather than absolute values.
    if omnibus.available
        omnibus.pvalueBoot = (1 + sum(Hboot >= omnibus.stat)) / ...
            (size(Hboot,1) + 1);
        omnibus.bootstrapAvailable = true;
        omnibus.bootstrapReason = '';
        omnibus.bootstrapCalibration = ...
            ['full-pipeline rank-selected ' bootstraptype ' bootstrap'];
        omnibus.statBoot = Hboot;
        omnibus.rankBoot = rankBoot;
        omnibus.relativeEigenvaluesBoot = relativeEigenvaluesBoot;
        omnibus.largestEigenvalueBoot = largestEigenvalueBoot;
        omnibus.rankThresholdBoot = rankThresholdBoot;
        omnibus.relativeRankThresholdBoot = relativeRankThresholdBoot;
        omnibus.rhoLastKeptBoot = rhoLastKeptBoot;
        omnibus.rhoFirstDiscardedBoot = rhoFirstDiscardedBoot;
        [rankValues,~,rankIndex] = unique(rankBoot);
        rankCounts = accumarray(rankIndex,1);
        rankFractions = rankCounts/numel(rankBoot);
        omnibus.rankBootFrequency = table(rankValues,rankCounts,rankFractions, ...
            'VariableNames',{'Rank','Count','Fraction'});
    end

    % Bootstrap confidence intervals.
    alphaCI = 1 - conflev;
    ciBoot = quantile(Tboot,[alphaCI/2 1-alphaCI/2],1);
end

% Human-readable table identifying explicitly the four monitored
% statistics and the analytical calibration used for each one. The bootstrap
% column is kept even when bootstrap=false so that the structure of
% out.results is invariant across calls.
statRowNames = {'TLmedian';'TLmean';'TDmedian';'TDmean'};
pvalueBootTable = NaN(4,1);
if bootstrap
    pvalueBootTable = pvalueBoot(:);
end

SEasym = NaN(4,1);
zAsym = NaN(4,1);
if isfield(asympt,'TLmean') && isstruct(asympt.TLmean)
    if isfield(asympt.TLmean,'se')
        SEasym(2) = asympt.TLmean.se;
    end
    if isfield(asympt.TLmean,'z')
        zAsym(2) = asympt.TLmean.z;
    end
end
if isfield(asympt,'TDmean') && isstruct(asympt.TDmean)
    if isfield(asympt.TDmean,'se')
        SEasym(4) = asympt.TDmean.se;
    end
    if isfield(asympt.TDmean,'z')
        zAsym(4) = asympt.TDmean.z;
    end
end
Calibration = local_calibration_labels(alpha,aggregation, ...
    consistencyfactor,coupledtrim);
results = table(Tobs(:),SEasym,zAsym,pvalueAsym(:),pvalueBootTable, ...
    Calibration,'VariableNames', ...
    {'Tobs','SEasym','zAsym','pvalueAsym','pvalueBoot','Calibration'}, ...
    'RowNames',statRowNames);

% Separate Gaussian-information benchmark. This is intentionally distinct
% from the primary alpha=0 ordinary calibration: TDmean is general-F and
% TLmean uses its empirical elliptical radial extension.
resultsGaussian = table();
if alpha == 0 && strcmp(aggregation,'ordinary') && ...
        isfield(asympt,'gaussian') && isstruct(asympt.gaussian)
    g = asympt.gaussian;
    gSE = [g.TLmean.se; g.TDmean.se];
    gz = [g.TLmean.z; g.TDmean.z];
    gp = [g.TLmean.pvalue; g.TDmean.pvalue];
    resultsGaussian = table([Tobs(2);Tobs(4)],gSE,gz,gp, ...
        'VariableNames',{'Tobs','SE','z','pvalue'}, ...
        'RowNames',{'TLmean';'TDmean'});
end

% Separate one-row table for the Hausman-type omnibus scatter discrepancy.
% The pvalue column retains the analytical chi-square calibration for backward
% compatibility; pvalueBoot is the full-pipeline rank-selected bootstrap value.
omnibusCalibration = {omnibus.calibration};
if omnibus.available
    resultsOmnibus = table(omnibus.stat,omnibus.df,omnibus.pvalue, ...
        omnibus.pvalueBoot,omnibusCalibration, ...
        'VariableNames',{'Stat','df','pvalue','pvalueBoot','Calibration'}, ...
        'RowNames',{'HausmanOmnibus'});
else
    resultsOmnibus = table(NaN,NaN,NaN,NaN,omnibusCalibration, ...
        'VariableNames',{'Stat','df','pvalue','pvalueBoot','Calibration'}, ...
        'RowNames',{'HausmanOmnibus'});
end

% Store output
out = struct;
out.pvalueAsym  = pvalueAsym;
out.Tobs        = Tobs;
out.results     = results;
out.resultsGaussian = resultsGaussian;
out.omnibus = omnibus;
out.resultsOmnibus = resultsOmnibus;
out.bootstrap   = bootstrap;
out.bootstraptype = bootstraptype;
out.dispresults = dispresults;
out.alpha       = alpha;
out.method      = method;
out.consistencyfactor = consistencyfactor;
out.adaptivepool = adaptivepool;
out.adaptiveminref = adaptiveminref;
out.omnibusrankexp = omnibusrankexp;
out.omnibusneff = omnibusneff;
if isfield(outFit,'kfactor')
    out.TEMkfactor = outFit.kfactor;
else
    out.TEMkfactor = NaN;
end
if isfield(outFit,'kinfo')
    out.TEMkinfo = outFit.kinfo;
else
    out.TEMkinfo = [];
end
out.aggregation = aggregation;
out.filterlev = filterlev;
out.filterCutoff = aggInfo.cutoff;
out.coupledtrim = coupledtrim;
out.coupledRows = find(forcedZeroTEM);
out.nCoupledCC = sum(forcedZeroTEM);
if isfield(outFit,'weights')
    out.TEMweights = outFit.weights(:);
    if isfield(outFit,'internalWeights')
        out.TEMinternalWeights = outFit.internalWeights(:);
    else
        out.TEMinternalWeights = outFit.weights(:);
    end
else
    % Untrimmed EM uses every row; expose this as unit TEM weights so that
    % the output fields have a consistent length for alpha=0 as well.
    out.TEMweights = true(n,1);
    out.TEMinternalWeights = true(n,1);
end
out.meanWeightsCC = meanWeightsCC;
completeRows = find(completeIdx);
out.meanRows = completeRows(meanWeightsCC);
out.nMeanCC = aggInfo.nSelected;
out.nComplete   = nComplete;
out.completeIdx = completeIdx;

if bootstrap
    out.pvalueBoot = pvalueBoot;
    out.Tboot = Tboot;
    out.ciBoot = ciBoot;
    out.nMeanBoot = nMeanBoot;
    out.nCoupledBoot = nCoupledBoot;
    out.nBootstrapAttempts = nBootstrapAttempts;
    out.nFailedBootstrapFits = nFailedBootstrapFits;
    out.nInvalidBootstrapStats = nInvalidBootstrapStats;
    out.nInvalidBootstrapOmnibus = nInvalidBootstrapOmnibus;
end

out.locCC = muCC;
out.covCC = SigCC;
out.robust = robust;
out.robustClass = robustClass;
out.robustEff = robustEff;
out.robustBonflev = robustBonflev;
out.outlierConflev = 0.99;

if alpha > 0 && strcmpi(robustClass,'FS')
    if isempty(robustBonflev)
        out.FSrule = 'resuperimposition';
    else
        out.FSrule = 'bonferroni';
    end
else
    out.FSrule = '';
end

if alpha > 0
    out.outliersCC = outRob.outliers;
    out.nOutliersCC = sum(local_cc_outlier_mask(outRob,nComplete));
else
    out.outliersCC = [];
    out.nOutliersCC = 0;
end

out.d2_cc       = d2_cc;
out.d2_all      = d2_all_cc;

out.loc       = muHat;
out.cov      = SigHat;
out.eps0 = eps0;
out.asympt = asympt;

% Optional commented display.
if dispresults
    switch aggregation
        case 'ordinary'
            procedure = 'MDP-U';
        case 'inlier'
            procedure = 'MDP-I';
        case 'fixed'
            procedure = 'MDP-F';
        otherwise
            procedure = 'MDP';
    end

    disp('*****************************************************************')
    disp('Mahalanobis Distance Perturbation test for MCAR')
    fprintf('Procedure                    : %s\n',procedure);
    fprintf('Mean aggregation             : %s\n',aggregation);
    fprintf('Trimming level alpha         : %g\n',alpha);
    fprintf('Complete cases               : %d out of %d\n',nComplete,n);
    fprintf('Complete cases used in means : %d\n',aggInfo.nSelected);
    if alpha == 0
        disp('Complete-case fit             : classical')
        disp('All-data fit                  : EM')
    else
        fprintf('Complete-case robust fit      : %s\n',robustClass);
        fprintf('All-data fit                  : TEM (%s consistency factor)\n', ...
            consistencyfactor);
    end
    fprintf('Distance rescaling method     : %s\n',method);
    if alpha == 0 && strcmp(aggregation,'ordinary')
        disp('Primary analytical reference  : statistic-specific; see Calibration column')
    elseif isfield(asympt,'mode') && ~isempty(asympt.mode)
        fprintf('Primary analytical reference  : %s\n',asympt.mode);
    end
    if bootstrap
        switch bootstraptype
            case 'gaussian'
                bootstrapLabel = 'Gaussian fitted complete-case model';
            case 'empiricalradial'
                bootstrapLabel = 'empirical-radial elliptical';
        end
        fprintf('Bootstrap null generator      : %s\n',bootstrapLabel);
        fprintf('Bootstrap replications kept   : %d\n',size(Tboot,1));
        fprintf('Bootstrap attempts            : %d\n',out.nBootstrapAttempts);
        if out.nFailedBootstrapFits > 0
            fprintf('Adaptive-TEM draws regenerated: %d\n', ...
                out.nFailedBootstrapFits);
        end
        if out.nInvalidBootstrapStats > 0
            fprintf('Nonfinite draws regenerated   : %d\n', ...
                out.nInvalidBootstrapStats);
        end
        if out.nInvalidBootstrapOmnibus > 0
            fprintf('Omnibus draws regenerated     : %d\n', ...
                out.nInvalidBootstrapOmnibus);
        end
    else
        disp('Bootstrap calibration          : not requested')
    end

    disp(' ')
    disp('Primary analytical results')
    disp(' ')
    disp(out.results)
    disp('Asymptotic inference is not available for the median statistics.')

    disp(' ')
    disp('Hausman-type omnibus scatter-discrepancy diagnostic')
    if out.omnibus.available
        disp(out.resultsOmnibus)
        fprintf('Estimated discrepancy rank    : %d of %d\n', ...
            out.omnibus.rank,out.omnibus.dimension);
        fprintf('Omnibus effective sample size : %g (%s)\n', ...
            out.omnibus.effectiveSampleSize, ...
            out.omnibus.effectiveSampleSizeSource);
        fprintf('Relative rank threshold eta_n : %.6g\n', ...
            out.omnibus.relativeRankThreshold);
        fprintf('Absolute eigenvalue threshold : %.6g\n', ...
            out.omnibus.rankThreshold);
        if isfinite(out.omnibus.rhoLastKept)
            fprintf('Last retained lambda/lambda1  : %.6g\n', ...
                out.omnibus.rhoLastKept);
        end
        if isfinite(out.omnibus.rhoFirstDiscarded)
            fprintf('First discarded lambda/lambda1: %.6g\n', ...
                out.omnibus.rhoFirstDiscarded);
        end
        if isfinite(out.omnibus.spectralGapRatio)
            fprintf('Boundary spectral gap ratio   : %.6g\n', ...
                out.omnibus.spectralGapRatio);
        end
        if isfinite(out.omnibus.lowerThresholdMargin)
            fprintf('Last-kept / eta margin        : %.6g\n', ...
                out.omnibus.lowerThresholdMargin);
        end
        if isfinite(out.omnibus.upperThresholdMargin)
            fprintf('Eta / first-discarded margin  : %.6g\n', ...
                out.omnibus.upperThresholdMargin);
        end
        if out.omnibus.nullSpaceResidual > 1e-6
            fprintf('Relative null-space residual  : %.6g\n', ...
                out.omnibus.nullSpaceResidual);
        end
        if bootstrap && out.omnibus.bootstrapAvailable
            fprintf('Omnibus bootstrap p-value     : %.6g\n', ...
                out.omnibus.pvalueBoot);
            disp('Bootstrap selected-rank frequencies')
            disp(out.omnibus.rankBootFrequency)
        elseif bootstrap && ~out.omnibus.bootstrapAvailable
            fprintf('Omnibus bootstrap unavailable : %s\n', ...
                out.omnibus.bootstrapReason);
        end
    else
        disp('Omnibus analytical calibration unavailable.')
        if ~isempty(out.omnibus.reason)
            fprintf('Reason: %s\n',out.omnibus.reason);
        end
    end

    if strcmp(aggregation,'inlier') && ...
            isfield(asympt,'inlierScaling') && ...
            isstruct(asympt.inlierScaling) && ...
            isfield(asympt.inlierScaling,'available') && ...
            asympt.inlierScaling.available
        fprintf('Stable-inlier retained fraction : %.6g\n', ...
            asympt.inlierScaling.fractionSelected);
        fprintf('MDP-I analytical direction      : generic retained moments\n');
        fprintf('TD retained first-moment norm   : %.6g\n', ...
            asympt.inlierScaling.firstMomentMahalanobis);
        fprintf('TD retained shape residual      : %.6g\n', ...
            asympt.inlierScaling.secondMomentShapeResidual);
        fprintf('TD radial special-case rI       : %.6g\n', ...
            asympt.inlierScaling.rI);
        fprintf('TL radial special-case rLI      : %.6g\n', ...
            asympt.inlierScaling.rLI);
    end

    % Surface analytical failures instead of leaving unexplained NaNs.
    meanAsymUnavailable = all(isnan(out.results.pvalueAsym([2 4])));
    if meanAsymUnavailable
        disp(' ')
        disp('WARNING: primary analytical calibration for the mean statistics is unavailable.')
        if isfield(asympt,'reason') && ~isempty(asympt.reason)
            fprintf('Reason: %s\n',asympt.reason);
        end
        if isfield(asympt,'rcondJA') && isfinite(asympt.rcondJA) && ...
                asympt.rcondJA <= 1e-12
            fprintf('rcond(J_A)                    : %.6g\n',asympt.rcondJA);
            disp(['The very small reciprocal condition number indicates severe ' ...
                'numerical ill-conditioning.'])
            disp(['This may be caused by variables on very different scales or ' ...
                'by near-linear dependence.'])
            disp('Consider rescaling the variables and recomputing the test.')
            disp('This numerical failure is not evidence for or against MCAR.')
        end
    end

    % For ordinary classical MDP-U, report the Gaussian-information benchmark
    % separately. It is not the primary calibration.
    if alpha == 0 && strcmp(aggregation,'ordinary')
        disp(' ')
        if isfield(asympt,'gaussian') && isstruct(asympt.gaussian) && ...
                isfield(asympt.gaussian,'available') && asympt.gaussian.available
            disp('Gaussian-information benchmark (normality assumed)')
            disp(' ')
            disp(out.resultsGaussian)
        elseif isfield(asympt,'gaussian') && isstruct(asympt.gaussian)
            disp('Gaussian-information benchmark unavailable.')
            if isfield(asympt.gaussian,'reason') && ~isempty(asympt.gaussian.reason)
                fprintf('Reason: %s\n',asympt.gaussian.reason);
            end
        end
    elseif alpha == 0 && strcmp(aggregation,'fixed')
        disp('Note: a general-F fixed-cutoff calibration is not currently available.')
    end

    if coupledtrim
        disp(['Note: analytical p-values use a frozen-eligibility approximation; ' ...
            'selection variability of the eligibility mask is not included.'])
    end
    disp('Small p-values indicate evidence against MCAR for the corresponding calibration.')
    disp('*****************************************************************')
end

% Optional plots
if plots
    statNames = {'median log-ratio', ...
        'mean log-ratio', ...
        'median difference', ...
        'mean difference'};

    figure;
    if bootstrap
        tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

        % Left part: scatter plot of the two distances.
        nexttile([2 1])
        scatter(d2_cc,d2_all_cc,'o')
        refline(1,0)
        xlabel('Complete-case distance')
        ylabel('Distance from all-data EM/TEM')
        title(['\alpha=' num2str(alpha) ', aggregation=' aggregation])
        box on

        % Right part: bootstrap distributions of the four statistics.
        for j = 1:4
            nexttile
            histogram(Tboot(:,j))
            hold on
            xline(Tobs(j),'r','LineWidth',1.5)
            title(['p_{boot}=' num2str(pvalueBoot(j),4)])
            xlabel(statNames{j})
            box on
        end
    else
        scatter(d2_cc,d2_all_cc,'o')
        refline(1,0)
        xlabel('Complete-case distance')
        ylabel('Distance from all-data EM/TEM')
        title(['\alpha=' num2str(alpha) ', aggregation=' aggregation])
        box on
    end
end

end

% -------------------------------------------------------------------------

function [muCC,SigCC,outRob] = ...
    local_complete_case_fit(Ycc,alpha,robustClass,robustEff,robustBonflev)
%local_complete_case_fit estimates location and scatter from complete cases.
%
% For all robust estimators, outliers are declared at simultaneous
% Bonferronized confidence level 0.99. For FSM, robustBonflev=0.99 uses
% the direct Bonferroni bound, whereas robustBonflev=[] uses the standard
% signal-detection, validation and envelope-resuperimposition rules. MCD
% and MM use the corresponding individual level 1-0.01/ncc.

ncc = size(Ycc,1);
conflevOut = 1-0.01/ncc;

if alpha == 0

    muCC = mean(Ycc,1);
    SigCC = cov(Ycc);
    outRob = [];

elseif strcmpi(robustClass,'FS')

    h = floor(ncc*(1-alpha));

    if isempty(robustBonflev)
        % Standard FSM signal detection and envelope resuperimposition.
        outRob = FSM(Ycc, ...
            'init',h, ...
            'plots',0, ...
            'msg',false);
    else
        % Direct Bonferroni-bound stopping rule.
        outRob = FSM(Ycc, ...
            'init',h, ...
            'bonflev',robustBonflev, ...
            'plots',0, ...
            'msg',false);
    end

    muCC = outRob.loc;
    SigCC = outRob.cov;

elseif strcmpi(robustClass,'MCD')

    outRob = mcd(Ycc, ...
        'bdp',alpha, ...
        'conflev',conflevOut, ...
        'plots',0, ...
        'msg',0);

    muCC = outRob.loc;
    SigCC = outRob.cov;

elseif strcmpi(robustClass,'MM')

    outRob = MMmult(Ycc, ...
        'Sbdp',alpha, ...
        'eff',robustEff, ...
        'conflev',conflevOut, ...
        'plots',0, ...
        'Smsg',0);

    muCC = outRob.loc;
    SigCC = outRob.cov;

else
    error('FSDA:mdMDPtest:WrongRobust', ...
        'Robust estimator must be ''FS'', ''MCD'' or ''MM''.');
end

end

% -------------------------------------------------------------------------
function [robustClass,robustEff,robustBonflev] = local_parse_robust(robust)
%local_parse_robust parses the robust option.

robustEff = 0.95;
robustBonflev = []; % default is envelope resuperimposition

if isstruct(robust)
    if ~isscalar(robust)
        error('FSDA:mdMDPtest:WrongRobust', ...
            'Option robust must be a scalar structure.');
    end

    allowedFields = {'class','eff','bonflev'};
    suppliedFields = fieldnames(robust);
    wrongFields = setdiff(suppliedFields,allowedFields);
    if ~isempty(wrongFields)
        error('FSDA:mdMDPtest:WrongRobustField', ...
            'Unknown field in option robust: %s.',wrongFields{1});
    end

    if ~isfield(robust,'class') || isempty(robust.class)
        error('FSDA:mdMDPtest:MissingRobustClass', ...
            ['When robust is a structure, field robust.class must be ' ...
            'specified as ''FS'', ''MCD'' or ''MM''.']);
    end
    robustClass = robust.class;

    if isfield(robust,'eff') && ~isempty(robust.eff)
        robustEff = robust.eff;
    end

    if isfield(robust,'bonflev')
        robustBonflev = robust.bonflev;
    end
else
    robustClass = robust;
end

if isstring(robustClass)
    if ~isscalar(robustClass)
        error('FSDA:mdMDPtest:WrongRobust', ...
            'Option robust.class must be a character vector or string scalar.');
    end
    robustClass = char(robustClass);
end

if ~ischar(robustClass) || size(robustClass,1) ~= 1
    error('FSDA:mdMDPtest:WrongRobust', ...
        'Option robust must specify ''FS'', ''MCD'' or ''MM''.');
end

robustClass = upper(strtrim(robustClass));
if ~any(strcmp(robustClass,{'FS','MCD','MM'}))
    error('FSDA:mdMDPtest:WrongRobust', ...
        'Option robust must specify ''FS'', ''MCD'' or ''MM''.');
end

if ~isscalar(robustEff) || ~isnumeric(robustEff) || ...
        ~isfinite(robustEff) || robustEff < 0.5 || robustEff > 0.99
    error('FSDA:mdMDPtest:WrongRobustEff', ...
        'Field robust.eff must be a scalar in the interval [0.5,0.99].');
end

if strcmp(robustClass,'FS')
    if ~(isempty(robustBonflev) || ...
            (isnumeric(robustBonflev) && isscalar(robustBonflev) && ...
            isfinite(robustBonflev) && robustBonflev == 0.99))
        error('FSDA:mdMDPtest:WrongRobustBonflev', ...
            ['Field robust.bonflev must be either 0.99 or empty. ' ...
            'Use 0.99 for the direct Bonferroni rule and [] for the ' ...
            'standard FSM envelope-resuperimposition rule.']);
    end
else
    if isstruct(robust) && isfield(robust,'bonflev')
        error('FSDA:mdMDPtest:WrongRobustBonflev', ...
            'Field robust.bonflev can be supplied only when robust.class=''FS''.');
    end
    robustBonflev = [];
end

end


function [d2_all_cc, muHat, SigHat, outFit] = local_fit_and_get_complete_distances( ...
    Y, completeIdx, alpha, method, tol, forcedZero, consistencyfactor, ...
    Yref, muref, sigmaref, adaptivepool, adaptiveminref)
%local_fit_and_get_complete_distances fits EM/TEM and returns complete-row distances.
%
% forcedZero is a logical n-vector. When alpha>0 and at least one entry is
% true, TEM's own concentration indicator is multiplied by ~forcedZero at
% every iteration. This implements the coupled sensitivity rule without
% changing the default mdTEM fit. For the adaptive correction, Yref, muref
% and sigmaref are the complete-case reference sample and fitted geometry.

p = size(Y,2);
n = size(Y,1);
if isempty(forcedZero)
    forcedZero = false(n,1);
else
    forcedZero = logical(forcedZero(:));
    if numel(forcedZero) ~= n
        error('FSDA:mdMDPtest:WrongCoupledMask', ...
            'The coupled TEM mask must have one entry for every row of Y.');
    end
end

if alpha == 0
    if any(forcedZero)
        error('FSDA:mdMDPtest:CoupledNeedsRobust', ...
            'Coupled TEM is available only when alpha>0.');
    end
    outFit = mdEM(Y);
elseif any(forcedZero)
    outFit = local_coupled_tem(Y,forcedZero,alpha,method,tol);
else
    temArgs = {'method',method,'alpha',alpha,'tol',tol, ...
        'consistencyfactor',consistencyfactor};
    if strcmp(consistencyfactor,'adaptive')
        temArgs = [temArgs, {'Yref',Yref,'muref',muref(:), ...
            'sigmaref',sigmaref,'adaptivepool',adaptivepool, ...
            'adaptiveminref',adaptiveminref}];
    end
    outFit = mdTEM(Y,temArgs{:});
    outFit.internalWeights = outFit.weights;
    outFit.forcedZero = forcedZero;
end

muHat = outFit.loc;
SigHat = outFit.cov;

[d2_part, poss] = mdPartialMD(Y, muHat, SigHat);
d2_full = mdPartialMD2full(d2_part, p, poss, 'method', method);
d2_all_cc = d2_full(completeIdx);

end

% -------------------------------------------------------------------------
function [nEff,source,spec] = local_resolve_omnibus_neff(spec,n,nComplete)
%local_resolve_omnibus_neff resolves the effective sample size for rank selection.
%
% The statistical rank threshold is eta_n=nEff^(-a). The choice of nEff
% affects only the spectral threshold used by the omnibus statistic.

if isstring(spec)
    if ~isscalar(spec)
        error('FSDA:mdMDPtest:WrongInputOpt', ...
            'Option omnibusneff must be ''n'', ''complete'' or a numeric scalar greater than one.');
    end
    spec = char(spec);
end

if ischar(spec)
    if size(spec,1) ~= 1
        error('FSDA:mdMDPtest:WrongInputOpt', ...
            'Option omnibusneff must be ''n'', ''complete'' or a numeric scalar greater than one.');
    end
    spec = lower(strtrim(spec));
    switch spec
        case 'n'
            nEff = n;
            source = 'n';
        case 'complete'
            nEff = nComplete;
            source = 'complete';
        otherwise
            error('FSDA:mdMDPtest:WrongInputOpt', ...
                'Option omnibusneff must be ''n'', ''complete'' or a numeric scalar greater than one.');
    end
elseif isnumeric(spec) && isscalar(spec) && isfinite(spec) && spec > 1
    nEff = double(spec);
    source = 'numeric';
else
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option omnibusneff must be ''n'', ''complete'' or a numeric scalar greater than one.');
end

if ~isfinite(nEff) || nEff <= 1
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'The effective sample size used by the omnibus rank rule must exceed one.');
end
end

% -------------------------------------------------------------------------
function consistencyfactor = local_parse_consistencyfactor(consistencyfactor)
%local_parse_consistencyfactor parses the mdTEM consistency-factor option.

if isstring(consistencyfactor)
    if ~isscalar(consistencyfactor)
        error('FSDA:mdMDPtest:WrongConsistencyFactor', ...
            ['Option consistencyfactor must be a character vector or ' ...
            'string scalar.']);
    end
    consistencyfactor = char(consistencyfactor);
end

if ~ischar(consistencyfactor) || size(consistencyfactor,1) ~= 1
    error('FSDA:mdMDPtest:WrongConsistencyFactor', ...
        ['Option consistencyfactor must be ''pattern'', ''adaptive'', ' ...
        '''global'', ''weighted'' or ''none''.']);
end

consistencyfactor = lower(strtrim(consistencyfactor));
valid = {'pattern','adaptive','global','weighted','none'};
if ~any(strcmp(consistencyfactor,valid))
    error('FSDA:mdMDPtest:WrongConsistencyFactor', ...
        ['Option consistencyfactor must be ''pattern'', ''adaptive'', ' ...
        '''global'', ''weighted'' or ''none''.']);
end
end

% -------------------------------------------------------------------------
function bootstraptype = local_parse_bootstraptype(bootstraptype)
%local_parse_bootstraptype parses the bootstrap null generator.

if isstring(bootstraptype)
    if ~isscalar(bootstraptype)
        error('FSDA:mdMDPtest:WrongBootstrapType', ...
            'Option bootstraptype must be a character vector or string scalar.');
    end
    bootstraptype = char(bootstraptype);
end

if ~ischar(bootstraptype) || size(bootstraptype,1) ~= 1
    error('FSDA:mdMDPtest:WrongBootstrapType', ...
        'Option bootstraptype must be ''gaussian'' or ''empiricalradial''.');
end

bootstraptype = lower(strtrim(bootstraptype));
if ~any(strcmp(bootstraptype,{'gaussian','empiricalradial'}))
    error('FSDA:mdMDPtest:WrongBootstrapType', ...
        'Option bootstraptype must be ''gaussian'' or ''empiricalradial''.');
end
end

% -------------------------------------------------------------------------
function aggregation = local_parse_aggregation(aggregation)
%local_parse_aggregation parses the mean aggregation option.

if isstring(aggregation)
    if ~isscalar(aggregation)
        error('FSDA:mdMDPtest:WrongAggregation', ...
            'Option aggregation must be a character vector or string scalar.');
    end
    aggregation = char(aggregation);
end


if ~ischar(aggregation) || size(aggregation,1) ~= 1
    error('FSDA:mdMDPtest:WrongAggregation', ...
        'Option aggregation must be ''auto'', ''ordinary'', ''inlier'' or ''fixed''.');
end


aggregation = lower(strtrim(aggregation));
if ~any(strcmp(aggregation,{'auto','ordinary','inlier','fixed'}))
    error('FSDA:mdMDPtest:WrongAggregation', ...
        'Option aggregation must be ''auto'', ''ordinary'', ''inlier'' or ''fixed''.');
end
end

% -------------------------------------------------------------------------
function [A,info] = local_aggregation_weights(d2_cc,outRob,alpha, ...
    aggregation,filterlev,p)
%local_aggregation_weights builds the complete-case mean aggregation rule.
%
% A is a logical vector indexed within the complete-case sample. The median
% statistics do not use A. For MDP-I, outlier labels are reconstructed from
% the complete-case robust fit. For MDP-F, the cutoff is fixed at the
% chi-square filterlev quantile and is applied to the complete-case distances.

ncc = numel(d2_cc);
A = true(ncc,1);
cutoff = NaN;

switch aggregation
    case 'ordinary'
        % MDP-U: all complete cases enter the two mean statistics.

    case 'inlier'
        % MDP-I: remove complete cases declared as outliers by the selected
        % robust complete-case estimator.
        if alpha == 0 || isempty(outRob) || ~isfield(outRob,'outliers')
            error('FSDA:mdMDPtest:MissingOutlierSet', ...
                ['Aggregation ''inlier'' requires a robust complete-case ' ...
                'fit returning the field outliers.']);
        end

        A = ~local_cc_outlier_mask(outRob,ncc);

    case 'fixed'
        % MDP-F: retain only complete cases inside the fixed robust-distance
        % region. The same fixed cutoff is used in every bootstrap sample.
        cutoff = chi2inv(filterlev,p);
        A = isfinite(d2_cc(:)) & d2_cc(:) <= cutoff;
end

info = struct;
info.nSelected = sum(A);
info.fractionSelected = info.nSelected/ncc;
info.cutoff = cutoff;
info.weights = A;
end

% -------------------------------------------------------------------------
function outlierMask = local_cc_outlier_mask(outRob,ncc)
%local_cc_outlier_mask converts robust complete-case outliers to a logical mask.

outlierMask = false(ncc,1);
if isempty(outRob) || ~isstruct(outRob) || ~isfield(outRob,'outliers')
    return
end

outliers = outRob.outliers;
if islogical(outliers) && numel(outliers) == ncc
    outlierMask = outliers(:);
elseif ~isempty(outliers)
    outliers = outliers(:);
    if isnumeric(outliers)
        outliers = outliers(isfinite(outliers) & outliers == floor(outliers) & ...
            outliers >= 1 & outliers <= ncc);
        outlierMask(unique(outliers)) = true;
    end
end
end

% -------------------------------------------------------------------------
function T = local_statistic(d2_cc, d2_all, eps0, meanWeights)
%local_statistic computes the four MDP test statistics.
%
% The median statistics use all complete cases. The two mean statistics use
% the logical selection vector meanWeights.

d2_cc = d2_cc(:);
d2_all = d2_all(:);
meanWeights = logical(meanWeights(:));

if numel(d2_cc) ~= numel(d2_all) || numel(meanWeights) ~= numel(d2_cc)
    error('FSDA:mdMDPtest:WrongAggregationWeights', ...
        'Distance vectors and aggregation weights must have the same length.');
end

rat = log((d2_all + eps0) ./ (d2_cc + eps0));
dif = d2_all - d2_cc;

T1 = median(rat);
T3 = median(dif);

if any(meanWeights)
    T2 = mean(rat(meanWeights));
    T4 = mean(dif(meanWeights));
else
    T2 = NaN;
    T4 = NaN;
end

T = [T1 T2 T3 T4];
end

% -------------------------------------------------------------------------
function [info,reason] = local_stable_inlier_scaling(Ycc,mu0,S0,weights,eps0)
%local_stable_inlier_scaling Generic retained-moment directions for MDP-I.
%
% Let X=Y-mu0, K0=S0^{-1}, gamma_I=P(omega=1), and define the retained
% moments
%
%   m_I = E(omega X)/gamma_I,
%   M_I = E(omega XX')/gamma_I.
%
% The generic first-order MDP-I directions are
%
%   b_I = 2*K0*m_I,
%   a_I = D_p'*vec(K0*M_I*K0).
%
% The stabilized log-ratio uses the analogous retained weighted moments
%
%   m_LI = E{omega X/(Q0+eps0)}/gamma_I,
%   M_LI = E{omega XX'/(Q0+eps0)}/gamma_I,
%
% with b_LI=2*K0*m_LI and a_LI=D_p'*vec(K0*M_LI*K0).
%
% The historical scalar radial multipliers rI and rLI are retained as
% diagnostics. They are the special proportional-moment case in which
% b_I=b_LI=0 and a_I=rI*a0, a_LI=rLI*a0.

info = struct('available',false,'genericAvailable',false,'reason','', ...
    'nComplete',NaN,'nSelected',NaN,'fractionSelected',NaN, ...
    'fractionExcluded',NaN,'rI',NaN,'rLI',NaN,'TLtoTDRatio',NaN, ...
    'meanQAll',NaN,'meanQInlier',NaN,'meanStabilizedQInlier',NaN, ...
    'kappa0',NaN,'kappaI',NaN,'lambdaI',NaN, ...
    'mI',[],'MI',[],'bI',[],'aI',[], ...
    'mLI',[],'MLI',[],'bLI',[],'aLI',[], ...
    'a0Radial',[],'aIRadial',[],'aLIRadial',[], ...
    'bINorm',NaN,'bLINorm',NaN, ...
    'aIRadialRelativeResidual',NaN,'aLIRadialRelativeResidual',NaN, ...
    'firstMomentMahalanobis',NaN,'logFirstMomentMahalanobis',NaN, ...
    'secondMomentShapeResidual',NaN,'logSecondMomentShapeResidual',NaN, ...
    'radialSigma2TD',NaN,'radialSigma2TL',NaN, ...
    'genericSigma2TD',NaN,'genericSigma2TL',NaN, ...
    'primaryCalibration','generic retained moments', ...
    'theoryStatus',['Generic retained-moment MDP-I calibration. The scalar ' ...
    'rI and rLI multipliers are retained only as the proportional-moment ' ...
    'special case.']);
reason = '';

if nargin == 0
    return
end

A = logical(weights(:));
[ncc,p] = size(Ycc);
info.nComplete = ncc;
info.nSelected = sum(A);
if ncc > 0
    info.fractionSelected = info.nSelected/ncc;
    info.fractionExcluded = 1-info.fractionSelected;
end

if numel(A) ~= ncc
    reason = ['The inlier aggregation weights and complete-case sample ' ...
        'have different lengths.'];
    return
end
if ncc == 0 || ~any(A)
    reason = ['Unable to evaluate stable-inlier retained moments because ' ...
        'no complete cases are retained.'];
    return
end
mu0 = mu0(:);
S0 = (S0+S0')/2;
if numel(mu0) ~= p || ~isequal(size(S0),[p p]) || ...
        any(~isfinite(mu0)) || any(~isfinite(S0(:))) || rcond(S0) <= 1e-12
    reason = 'The complete-case reference geometry is invalid or singular.';
    return
end
if ~isscalar(eps0) || ~isfinite(eps0) || eps0 <= 0
    reason = 'The log-ratio stabilizer eps0 must be finite and positive.';
    return
end

K0 = S0\eye(p);
K0 = (K0+K0')/2;
X = Ycc-mu0';
q = sum((X*K0).*X,2);
if any(~isfinite(q)) || any(q < 0)
    reason = ['Unable to evaluate stable-inlier retained moments because ' ...
        'the complete-case squared distances are negative or nonfinite.'];
    return
end

XA = X(A,:);
qA = q(A);
nA = size(XA,1);
meanQAll = mean(q);
meanQInlier = mean(qA);
meanStabilizedQInlier = mean(qA./(qA+eps0));
if ~isfinite(meanQAll) || meanQAll <= 0 || ...
        ~isfinite(meanQInlier) || meanQInlier <= 0 || ...
        ~isfinite(meanStabilizedQInlier) || meanStabilizedQInlier <= 0
    reason = ['The empirical stable-inlier radial moments are not finite ' ...
        'and positive.'];
    return
end

mI = mean(XA,1)';
MI = (XA'*XA)/nA;
MI = (MI+MI')/2;
den = qA+eps0;
mLI = mean(bsxfun(@rdivide,XA,den),1)';
MLI = XA'*bsxfun(@rdivide,XA,den)/nA;
MLI = (MLI+MLI')/2;

Dp = local_duplication_matrix(p);
KMK = K0*MI*K0;
KMLK = K0*MLI*K0;
bI = 2*K0*mI;
aI = Dp'*KMK(:);
bLI = 2*K0*mLI;
aLI = Dp'*KMLK(:);

kappa0 = meanQAll/p;
kappaI = meanQInlier/p;
lambdaI = meanStabilizedQInlier/p;
rI = kappaI/kappa0;
rLI = lambdaI/kappa0;
a0Radial = kappa0*(Dp'*K0(:));
aIRadial = rI*a0Radial;
aLIRadial = rLI*a0Radial;

R0 = chol(S0,'lower');
MIstd = R0\MI/R0';
MLIstd = R0\MLI/R0';
shapeI = MIstd-(trace(MIstd)/p)*eye(p);
shapeLI = MLIstd-(trace(MLIstd)/p)*eye(p);
firstMah = sqrt(max(real(mI'*K0*mI),0));
logFirstMah = sqrt(max(real(mLI'*K0*mLI),0));

info.available = true;
info.genericAvailable = true;
info.rI = rI;
info.rLI = rLI;
info.TLtoTDRatio = rLI/rI;
info.meanQAll = meanQAll;
info.meanQInlier = meanQInlier;
info.meanStabilizedQInlier = meanStabilizedQInlier;
info.kappa0 = kappa0;
info.kappaI = kappaI;
info.lambdaI = lambdaI;
info.mI = mI;
info.MI = MI;
info.bI = bI;
info.aI = aI;
info.mLI = mLI;
info.MLI = MLI;
info.bLI = bLI;
info.aLI = aLI;
info.a0Radial = a0Radial;
info.aIRadial = aIRadial;
info.aLIRadial = aLIRadial;
info.bINorm = norm(bI);
info.bLINorm = norm(bLI);
info.aIRadialRelativeResidual = norm(aI-aIRadial)/max(norm(aI),eps);
info.aLIRadialRelativeResidual = norm(aLI-aLIRadial)/max(norm(aLI),eps);
info.firstMomentMahalanobis = firstMah;
info.logFirstMomentMahalanobis = logFirstMah;
info.secondMomentShapeResidual = norm(shapeI,'fro')/max(norm(MIstd,'fro'),eps);
info.logSecondMomentShapeResidual = norm(shapeLI,'fro')/max(norm(MLIstd,'fro'),eps);
end

% -------------------------------------------------------------------------
function asympt = local_apply_aggregation_asymptotic(asympt,aggregation, ...
    aggInfo,p,eps0,n,Tobs,d2cc,Ycc,muCC,SigCC)
%local_apply_aggregation_asymptotic applies MDP-U/I/F first-order calibration.
%
% MDP-U and MDP-F retain scalar first-order coefficients. MDP-I instead uses
% the generic stable-inlier directions bI,aI and bLI,aLI together with the
% joint location--scatter estimator-discrepancy covariance. The empirical
% radial ratios rI and rLI are preserved as the proportional-moment special
% case. For MDP-F, c is fixed as n grows.

asympt.aggregation = aggregation;
asympt.filterCutoff = NaN;
asympt.kc = NaN;
asympt.lambda = NaN;
asympt.gammaFilter = NaN;
[asympt.inlierScaling,~] = local_stable_inlier_scaling();

% The stable-inlier multipliers can be evaluated even if the base sandwich
% is unavailable. This preserves useful diagnostics when analytical
% studentization fails for another reason.
if strcmp(aggregation,'inlier')
    if ~isfield(aggInfo,'weights')
        asympt.available = false;
        asympt.reason = ['The estimator-induced inlier indicators are not ' ...
            'available for analytical scaling.'];
        return
    end
    [asympt.inlierScaling,scaleReason] = ...
        local_stable_inlier_scaling(Ycc,muCC,SigCC,aggInfo.weights,eps0);
    asympt.inlierScaling.reason = scaleReason;
    asympt.kc = asympt.inlierScaling.rI;
    asympt.lambda = asympt.inlierScaling.rLI;
    asympt.gammaFilter = asympt.inlierScaling.fractionSelected;
    if ~isempty(scaleReason)
        asympt.available = false;
        asympt.reason = scaleReason;
        return
    end
end

% If the base estimator-discrepancy benchmark is unavailable, retain its
% reason and only document the requested aggregation.
if ~isfield(asympt,'available') || ~asympt.available
    if strcmp(aggregation,'fixed')
        asympt.filterCutoff = aggInfo.cutoff;
    end
    return
end

baseSigmaD2 = asympt.TDmean.sigma2;
baseDegenerate = asympt.degenerate;
if ~isfinite(baseSigmaD2) || baseSigmaD2 < 0
    asympt.available = false;
    asympt.reason = 'The base analytical variance is not finite.';
    return
end

% Store the unscaled estimator-discrepancy variance for diagnostics.
asympt.baseSigmaD2 = baseSigmaD2;
if isfield(asympt,'TLmean') && isfield(asympt.TLmean,'kappa')
    baseKappa = asympt.TLmean.kappa;
else
    baseKappa = NaN;
end
asympt.baseKappa = baseKappa;

% Generic stable-inlier calibration. The primary MDP-I law uses the joint
% location--scatter estimator-discrepancy covariance and the retained-moment
% directions (bI,aI) and (bLI,aLI). The scalar rI/rLI laws are retained only
% as proportional-moment diagnostics.
if strcmp(aggregation,'inlier')
    sdim = p*(p+1)/2;
    if ~isfield(asympt,'OmegaJoint') || ~isnumeric(asympt.OmegaJoint) || ...
            ~isequal(size(asympt.OmegaJoint),[p+sdim p+sdim]) || ...
            any(~isfinite(asympt.OmegaJoint(:)))
        asympt.available = false;
        if isfield(asympt,'jointReason') && ~isempty(asympt.jointReason)
            asympt.reason = asympt.jointReason;
        else
            asympt.reason = ['The joint location--scatter estimator-discrepancy ' ...
                'covariance required for generic MDP-I calibration is unavailable.'];
        end
        return
    end

    dirD = [asympt.inlierScaling.bI; asympt.inlierScaling.aI];
    dirL = [asympt.inlierScaling.bLI; asympt.inlierScaling.aLI];
    OmegaJoint = (asympt.OmegaJoint+asympt.OmegaJoint')/2;
    sigmaD2Generic = real(dirD'*OmegaJoint*dirD);
    sigmaL2Generic = real(dirL'*OmegaJoint*dirL);

    omegaScale = norm(OmegaJoint,2);
    tolD = 1e-10*max(1,omegaScale*max(1,norm(dirD)^2));
    tolL = 1e-10*max(1,omegaScale*max(1,norm(dirL)^2));
    if sigmaD2Generic < 0 && sigmaD2Generic >= -tolD
        sigmaD2Generic = 0;
    end
    if sigmaL2Generic < 0 && sigmaL2Generic >= -tolL
        sigmaL2Generic = 0;
    end
    if ~isfinite(sigmaD2Generic) || sigmaD2Generic < 0 || ...
            ~isfinite(sigmaL2Generic) || sigmaL2Generic < 0
        asympt.available = false;
        asympt.reason = ['The generic retained-moment MDP-I variance is not ' ...
            'finite and nonnegative.'];
        return
    end

    % Preserve the old scalar radial calibration as an explicit special-case
    % benchmark, but do not use it for the primary p-values.
    radialSigma2TD = asympt.inlierScaling.rI^2*baseSigmaD2;
    radialSigma2TL = asympt.inlierScaling.rLI^2*baseSigmaD2;
    asympt.inlierScaling.radialSigma2TD = radialSigma2TD;
    asympt.inlierScaling.radialSigma2TL = radialSigma2TL;
    asympt.inlierScaling.genericSigma2TD = sigmaD2Generic;
    asympt.inlierScaling.genericSigma2TL = sigmaL2Generic;

    if isfield(asympt.TDmean,'components') && isstruct(asympt.TDmean.components)
        asympt.TDmean.baseComponents = asympt.TDmean.components;
        asympt.TDmean = rmfield(asympt.TDmean,'components');
    end
    asympt.TDmean.coefficient = NaN;
    asympt.TDmean.coefficientType = 'generic retained-moment joint direction';
    asympt.TDmean.bI = asympt.inlierScaling.bI;
    asympt.TDmean.aI = asympt.inlierScaling.aI;
    asympt.TDmean.rI = asympt.inlierScaling.rI;
    asympt.TDmean.radialSpecialCaseSigma2 = radialSigma2TD;
    asympt.TDmean.sigma2 = sigmaD2Generic;
    asympt.TDmean.se = sqrt(sigmaD2Generic/n);

    asympt.TLmean.coefficient = NaN;
    asympt.TLmean.coefficientType = 'generic retained weighted-moment joint direction';
    asympt.TLmean.bLI = asympt.inlierScaling.bLI;
    asympt.TLmean.aLI = asympt.inlierScaling.aLI;
    asympt.TLmean.rLI = asympt.inlierScaling.rLI;
    asympt.TLmean.TLtoTDRatio = asympt.inlierScaling.TLtoTDRatio;
    asympt.TLmean.radialSpecialCaseSigma2 = radialSigma2TL;
    if ~isfield(asympt.TLmean,'kappa')
        asympt.TLmean.kappa = baseKappa;
    end
    asympt.TLmean.lambda = asympt.inlierScaling.rLI;
    asympt.TLmean.sigma2 = sigmaL2Generic;
    asympt.TLmean.se = sqrt(sigmaL2Generic/n);

    % Backward-compatible scalar fields now document the radial special case.
    asympt.kc = asympt.inlierScaling.rI;
    asympt.lambda = asympt.inlierScaling.rLI;
    asympt.gammaFilter = asympt.inlierScaling.fractionSelected;

    asympt.TDmean.degenerate = sigmaD2Generic <= tolD;
    asympt.TLmean.degenerate = sigmaL2Generic <= tolL;
    asympt.degenerate = asympt.TDmean.degenerate;

    if asympt.TDmean.degenerate
        asympt.TDmean.z = NaN;
        asympt.TDmean.pvalue = NaN;
    else
        asympt.TDmean.z = Tobs(4)/asympt.TDmean.se;
        asympt.TDmean.pvalue = erfc(abs(asympt.TDmean.z)/sqrt(2));
    end
    if asympt.TLmean.degenerate
        asympt.TLmean.z = NaN;
        asympt.TLmean.pvalue = NaN;
    else
        asympt.TLmean.z = Tobs(2)/asympt.TLmean.se;
        asympt.TLmean.pvalue = erfc(abs(asympt.TLmean.z)/sqrt(2));
    end

    note = ['Aggregation=''inlier'' uses the generic retained-moment ' ...
        'joint location--scatter calibration. The empirical rI and rLI ' ...
        'multipliers are reported only as the retained-moment proportionality ' ...
        'special case. A nonvanishing excluded fraction is permitted under ' ...
        'the stable-inlier conditions.'];
    if isfield(asympt,'theoryStatus') && ~isempty(asympt.theoryStatus)
        asympt.theoryStatus = [asympt.theoryStatus ' ' note];
    else
        asympt.theoryStatus = note;
    end
    return
end

switch aggregation
    case 'ordinary'
        coefD = 1;
        coefL = baseKappa;
        asympt.kc = 1;
        asympt.lambda = baseKappa;
        asympt.gammaFilter = 1;

    case 'fixed'
        c = aggInfo.cutoff;
        isAdaptiveRadial = isfield(asympt,'radialModel') && ...
            strcmp(asympt.radialModel,'adaptive-empirical');

        if isAdaptiveRadial
            qref = d2cc(:);
            qref = qref(isfinite(qref) & qref >= 0);
            inside = qref <= c;
            gamma = mean(inside);
            meanQ = mean(qref);
            if isempty(qref) || ~isfinite(c) || c <= 0 || ...
                    ~isfinite(gamma) || gamma <= 0 || ...
                    ~isfinite(meanQ) || meanQ <= 0 || ~any(inside)
                asympt.available = false;
                asympt.reason = ['Unable to evaluate the adaptive fixed-cutoff ' ...
                    'radial coefficients from the complete-case distances.'];
                return
            end

            % Under the elliptical common-target model S0=tau*Sigma,
            % the fixed-cutoff TD coefficient relative to MDP-U is
            % E(Q0|Q0<=c)/E(Q0). The stabilized log-ratio coefficient is
            % E{Q0/(Q0+eps0)|Q0<=c}/E(Q0).
            kc = mean(qref(inside))/meanQ;
            lambda = mean(qref(inside)./(qref(inside)+eps0))/meanQ;
            note = ['Aggregation=''fixed'' uses empirical elliptical radial ' ...
                'coefficients estimated from the robust complete-case ' ...
                'distances.'];
        else
            gamma = chi2cdf(c,p);
            if ~isfinite(c) || c <= 0 || ~isfinite(gamma) || gamma <= 0
                asympt.available = false;
                asympt.reason = ['Unable to evaluate the fixed-cutoff analytical ' ...
                    'coefficients.'];
                return
            end
            kc = chi2cdf(c,p+2)/gamma;
            lambda = integral(@(q) local_kappa_integrand(q,p,eps0),0,c, ...
                'RelTol',1e-10,'AbsTol',1e-12)/(p*gamma);
            note = ['Aggregation=''fixed'' applies the fixed-cutoff first-order ' ...
                'Tallis scaling to the base estimator-discrepancy variance.'];
        end

        if ~isfinite(kc) || kc <= 0 || ~isfinite(lambda) || lambda <= 0
            asympt.available = false;
            asympt.reason = ['The fixed-cutoff analytical coefficients are ' ...
                'not finite and positive.'];
            return
        end

        coefD = kc;
        coefL = lambda;
        asympt.filterCutoff = c;
        asympt.kc = kc;
        asympt.lambda = lambda;
        asympt.gammaFilter = gamma;

        if isfield(asympt,'theoryStatus') && ~isempty(asympt.theoryStatus)
            asympt.theoryStatus = [asympt.theoryStatus ' ' note];
        else
            asympt.theoryStatus = note;
        end

    otherwise
        error('FSDA:mdMDPtest:WrongAggregation', ...
            'Unknown aggregation rule in analytical scaling.');
end

if ~isfinite(coefD) || coefD <= 0
    asympt.available = false;
    asympt.reason = 'The analytical distance-difference coefficient is not finite.';
    return
end
if ~isfinite(coefL) || coefL <= 0
    asympt.available = false;
    asympt.reason = 'The analytical log-ratio coefficient is not finite.';
    return
end

% Apply the selected scalar coefficient to the common estimator-discrepancy
% influence law.
asympt.TDmean.coefficient = coefD;
asympt.TDmean.sigma2 = coefD^2*baseSigmaD2;
asympt.TDmean.se = sqrt(asympt.TDmean.sigma2/n);

asympt.TLmean.coefficient = coefL;
asympt.TLmean.sigma2 = coefL^2*baseSigmaD2;
asympt.TLmean.se = sqrt(asympt.TLmean.sigma2/n);

% kappa remains the ordinary stabilized-log coefficient. The field lambda
% is the aggregation-specific coefficient actually used.
if ~isfield(asympt.TLmean,'kappa')
    asympt.TLmean.kappa = baseKappa;
end
asympt.TLmean.lambda = coefL;
if strcmp(aggregation,'inlier')
    asympt.TDmean.rI = coefD;
    asympt.TLmean.rLI = coefL;
    asympt.TLmean.TLtoTDRatio = coefL/coefD;
end

% TEM variance components, when available, inherit the same TD multiplier.
if isfield(asympt.TDmean,'components') && isstruct(asympt.TDmean.components)
    fn = fieldnames(asympt.TDmean.components);
    for j = 1:numel(fn)
        val = asympt.TDmean.components.(fn{j});
        if isnumeric(val) && isscalar(val) && isfinite(val)
            asympt.TDmean.components.(fn{j}) = coefD^2*val;
        end
    end
end

% Preserve the base TEM and complete-case diagnostics for backward
% compatibility, and add their contribution on the selected TD scale.
if isfield(asympt,'TEM') && isstruct(asympt.TEM) && ...
        isfield(asympt.TEM,'sigma2') && isfinite(asympt.TEM.sigma2)
    asympt.TEM.sigma2Selected = coefD^2*asympt.TEM.sigma2;
end
if isfield(asympt,'completeCase') && isstruct(asympt.completeCase)
    if isfield(asympt.completeCase,'sigma2') && ...
            isfinite(asympt.completeCase.sigma2)
        asympt.completeCase.sigma2Selected = ...
            coefD^2*asympt.completeCase.sigma2;
    end
    if isfield(asympt.completeCase,'crossTEMCC') && ...
            isfinite(asympt.completeCase.crossTEMCC)
        asympt.completeCase.crossTEMCCSelected = ...
            coefD^2*asympt.completeCase.crossTEMCC;
    end
end

% Degeneracy is unchanged by a strictly positive scalar multiplier. Keep
% the numerical decision made by the base analytical routine.
asympt.degenerate = baseDegenerate;

if asympt.degenerate
    asympt.TDmean.z = NaN;
    asympt.TDmean.pvalue = NaN;
    asympt.TLmean.z = NaN;
    asympt.TLmean.pvalue = NaN;
else
    asympt.TDmean.z = Tobs(4)/asympt.TDmean.se;
    asympt.TDmean.pvalue = erfc(abs(asympt.TDmean.z)/sqrt(2));
    asympt.TLmean.z = Tobs(2)/asympt.TLmean.se;
    asympt.TLmean.pvalue = erfc(abs(asympt.TLmean.z)/sqrt(2));
end
end

% -------------------------------------------------------------------------
function asympt = local_tem_adaptive_asymptotic(Y, maskMiss, completeIdx, ...
    outFit, alpha, method, robustClass, robustEff, robustBonflev, ...
    Tobs, eps0, muCC, SigCC, adaptivepool, adaptiveminref, outRob)
%local_tem_adaptive_asymptotic Adaptive alpha>0 TEM sandwich calibration.
%
% This routine implements the explicit first-order expansion for the
% distribution-adaptive radial consistency correction under the local shell-
% moment conditions of the paper. Ellipticity is sufficient, but not necessary,
% for these identities. The complete-case robust scatter converges to a common
% target S0, not necessarily to the covariance matrix itself.
%
% With adaptivepool=false, pattern g uses
%
%   kappa_g = mean(Q_g | Q_g<=a_g)/p_g.
%
% With adaptivepool=true, which is the default, projected complete-case
% radii are pooled across patterns having the same observed dimension k.
% The pooled factor is
%
%                 sum_i C_i sum_{g in G_k} Q_ig I(Q_ig<=a_k)
%   kappa_k = ------------------------------------------------------- .
%              k sum_i C_i sum_{g in G_k} I(Q_ig<=a_k)
%
% The several Q_ig generated by one complete row are generally dependent.
% Consequently the pooled factor influence is constructed at the complete-
% observation level; the m_k*ncc projections are never treated as independent
% reference observations. The random common TEM cutoff is absorbed through
% the effective Jacobian in both the pooled and unpooled cases.
%
% The boundary density is estimated by a Gaussian kernel density estimator
% with reflection at zero. The analytical Jacobian covers the parameter-free
% adjusted-distance mappings 'pri', 'expScale', 'zMap', 'chiMap' and
% 'betaMap'. Each mapping enters through its inverse raw-distance cutoff
% a_g(c). For all five mappings, patterns having the same observed dimension
% have the same raw cutoff, as required by dimension pooling. The
% scatter-dependent mapping 'detMap' is deliberately excluded.
%
% Although the common empirical TEM cutoff c is random, no derivative of the
% map with respect to c is needed in the first-order scatter Jacobian. At the
% common target, the consistency-corrected population estimating equation is
% zero for every admissible raw cutoff a_g. Hence the direct moving-boundary
% term and the cutoff-induced change of the consistency factor cancel.

[n,p] = size(Y);
qhat = sum(completeIdx)/n;
s = p*(p+1)/2;

asympt = local_unavailable_asymptotic(qhat);
asympt.reason = '';
asympt.VA = NaN;
if adaptivepool
    asympt.mode = 'dimension-pooled adaptive TEM influence-function sandwich';
else
    asympt.mode = 'adaptive TEM influence-function sandwich';
end
asympt.radialModel = 'adaptive-empirical';
asympt.theoryStatus = ['available under the local shell-moment conditions ' ...
    'for the explicit adaptive TEM sandwich; ellipticity is sufficient but not necessary'];
asympt.TEM = struct('available',false,'sigma2',NaN,'threshold',NaN, ...
    'Jrcond',NaN,'Jmurcond',NaN,'nActivePatterns',0,'nPatterns',0, ...
    'factorIFIncluded',true,'densityMethod','Gaussian KDE with reflection', ...
    'adaptivepool',adaptivepool,'nFactorGroups',0, ...
    'sigma2Base',NaN,'sigma2Factor',NaN,'twiceBaseFactor',NaN, ...
    'meanInfluenceBeforeCentering',NaN,'meanBaseBeforeCentering',NaN, ...
    'meanFactorBeforeCentering',NaN,'estimatingEquationMeanNorm',NaN, ...
    'patternDiagnostics',table(),'factorDiagnostics',table());
asympt.completeCase = struct('influenceMethod','', ...
    'sigma2',NaN,'crossTEMCC',NaN,'nComplete',sum(completeIdx), ...
    'nReferenceInliers',NaN,'nReferenceOutliers',NaN, ...
    'frozenClassification',false,'robustBonflev',robustBonflev, ...
    'relCovFrozenVsFS',NaN,'relLocFrozenVsFS',NaN, ...
    'meanInfluenceBeforeCentering',NaN);

supportedMethods = {'pri','expScale','zMap','chiMap','betaMap'};
if strcmpi(method,'detMap')
    asympt.reason = ['For consistencyfactor=''adaptive'' the analytical TEM ' ...
        'sandwich deliberately excludes method=''detMap'' because its ' ...
        'adjusted-distance map depends explicitly on the scatter matrix.'];
    return
elseif ~any(strcmpi(method,supportedMethods))
    asympt.reason = ['For consistencyfactor=''adaptive'' the analytical TEM ' ...
        'sandwich is implemented for ''pri'', ''expScale'', ''zMap'', ' ...
        '''chiMap'' and ''betaMap''.'];
    return
end

if isempty(outFit) || ~isfield(outFit,'weights') || ...
        ~isfield(outFit,'loc') || ~isfield(outFit,'cov') || ...
        ~isfield(outFit,'kinfo')
    asympt.reason = 'The adaptive TEM object does not contain the required fields.';
    return
end

Ycc = Y(completeIdx,:);
ncc = size(Ycc,1);

% The adaptive radial factors are estimated in the frozen complete-case
% reference geometry. Keep this geometry separate from the fitted TEM
% geometry used to evaluate the TEM estimating equation below. Both are
% consistent for the same population target under the null, but using the
% fitted TEM geometry for the TEM-side innovation makes the empirical
% estimating equation self-consistent in finite samples.
muRef = muCC(:);
Sref = (SigCC+SigCC')/2;
if numel(muRef) ~= p || ~isequal(size(Sref),[p p]) || ...
        any(~isfinite(Sref(:))) || rcond(Sref) <= 1e-12
    asympt.reason = 'The robust complete-case reference geometry is singular.';
    return
end

muTEM = outFit.loc(:);
STEM = (outFit.cov+outFit.cov')/2;
if numel(muTEM) ~= p || ~isequal(size(STEM),[p p]) || ...
        any(~isfinite(STEM(:))) || rcond(STEM) <= 1e-12
    asympt.reason = 'The fitted adaptive TEM geometry is singular.';
    return
end

% The ordinary MDP-U coefficient at a common target S0=tau*Sigma is
% a0=(1/tau)D_p' vec(S0^{-1}). Estimate 1/tau by mean(Q0)/p, where Q0 is
% the robust complete-case squared distance. This avoids assuming S0=Sigma.
Dp = local_duplication_matrix(p);
K0 = Sref\eye(p);
K0 = (K0+K0')/2;
q0 = mahalFS(Ycc,muRef',Sref);
q0 = q0(:);
if any(~isfinite(q0)) || any(q0 < 0)
    asympt.reason = 'Nonfinite robust complete-case distances in adaptive calibration.';
    return
end
meanQ0 = mean(q0);
if ~isfinite(meanQ0) || meanQ0 <= 0
    asympt.reason = 'The mean robust complete-case squared distance is not positive.';
    return
end
invTauHat = meanQ0/p;
a0 = invTauHat*(Dp'*K0(:));

% Under ellipticity, the stabilized mean log-ratio coefficient is a scalar
% multiple of the MDP-U coefficient. On the S0 radial scale this multiplier
% is E{Q0/(Q0+eps0)}/E(Q0).
baseKappa = mean(q0./(q0+eps0))/meanQ0;
if ~isfinite(baseKappa) || baseKappa <= 0
    asympt.reason = 'The adaptive stabilized-log coefficient is not finite.';
    return
end

% Reconstruct the common adjusted-distance threshold. mdTEM stores the
% final threshold explicitly. The adjusted-distance fallback is retained for
% compatibility with fitted objects created by older versions.
cthr = NaN;
if isfield(outFit,'cthr') && isscalar(outFit.cthr) && isfinite(outFit.cthr)
    cthr = outFit.cthr;
end

if isfield(outFit,'adjustedD2') && numel(outFit.adjustedD2)==n
    d2adj = outFit.adjustedD2(:);
else
    d2adj = local_adjusted_distances(Y,outFit.loc,outFit.cov,method);
end
w = outFit.weights(:) > 0;
if ~isfinite(cthr)
    if ~any(w & isfinite(d2adj))
        asympt.reason = 'Unable to reconstruct the adaptive TEM threshold.';
        return
    end
    cthr = max(d2adj(w & isfinite(d2adj)));
end
asympt.TEM.threshold = cthr;

% Full complete-case scatter influence is needed because each empirical
% radial factor depends on the robust reference scatter through Q_g. For FS
% this is evaluated analytically conditional on the full-sample final FS
% classification. MCD and MM retain the delete-one jackknife for now.
[PsiMuCcc,PsiCcc,ccInfo,ccReason] = local_cc_joint_influence(Ycc,alpha, ...
    robustClass,robustEff,robustBonflev,outRob);
if ~isempty(ccReason)
    asympt.reason = ccReason;
    return
end
if size(PsiMuCcc,2) ~= p || size(PsiCcc,2) ~= s
    asympt.reason = ['The robust complete-case location/scatter influence has ' ...
        'the wrong dimension.'];
    return
end
asympt.completeCase.influenceMethod = ccInfo.influenceMethod;
asympt.completeCase.nReferenceInliers = ccInfo.nReferenceInliers;
asympt.completeCase.nReferenceOutliers = ccInfo.nReferenceOutliers;
asympt.completeCase.frozenClassification = ccInfo.frozenClassification;
asympt.completeCase.relCovFrozenVsFS = ccInfo.relCovFrozenVsFS;
asympt.completeCase.relLocFrozenVsFS = ccInfo.relLocFrozenVsFS;

PsiCfull = zeros(n,s);
PsiCfull(completeIdx,:) = PsiCcc/qhat;
PsiMuCfull = zeros(n,p);
PsiMuCfull(completeIdx,:) = PsiMuCcc/qhat;

% Missingness patterns and empirical pattern probabilities.
[patt,~,ic] = unique(maskMiss,'rows');
G = size(patt,1);
asympt.TEM.nPatterns = G;

J = zeros(s,s);
Jmu = zeros(p,p);
XiBase = zeros(n,s);
XiMu = zeros(n,p);
XiFactor = zeros(n,s);
active = false(G,1);

pobsDiag = zeros(G,1);
piDiag = zeros(G,1);
nkeptDiag = zeros(G,1);
aDiag = NaN(G,1);
HDiag = NaN(G,1);
kDiag = NaN(G,1);
fDiag = NaN(G,1);
bwDiag = NaN(G,1);
rhoDiag = NaN(G,1);
nrefDiag = zeros(G,1);
nrefkeptDiag = zeros(G,1);
poolSizeDiag = zeros(G,1);
factorGroupDiag = zeros(G,1);

obsCell = cell(G,1);
Lcell = cell(G,1);
Vcell = cell(G,1);
Mcell = cell(G,1);
BgTEMcell = cell(G,1);
QrefCell = cell(G,1);
bgMat = NaN(s,G);
kVec = NaN(G,1);

% First pass: construct pattern-specific reference radii and TEM geometry.
% Radial factors are computed in a second pass so that equal-dimensional
% patterns can share one pooled factor when adaptivepool=true.
for g = 1:G
    obs = find(~patt(g,:));
    mis = find(patt(g,:));
    pg = numel(obs);
    pig = sum(ic==g)/n;
    nkeptg = sum(w & ic==g);

    pobsDiag(g) = pg;
    piDiag(g) = pig;
    nkeptDiag(g) = nkeptg;
    obsCell{g} = obs;

    if pg == 0
        continue
    end

    ag = local_tem_inv_adjust(cthr,pg,p,n,method);
    aDiag(g) = ag;
    if ~isfinite(ag) || ag <= 0 || nkeptg == 0
        continue
    end

    % Frozen complete-case reference geometry.
    SgRef = Sref(obs,obs);
    SgRef = (SgRef+SgRef')/2;
    if rcond(SgRef) <= 1e-12
        asympt.reason = ['A reference pattern scatter block is numerically ' ...
            'singular in the adaptive calibration.'];
        return
    end
    BgRef = SgRef\eye(pg);
    BgRef = (BgRef+BgRef')/2;

    Zref = Ycc(:,obs)-muRef(obs)';
    QrefCell{g} = sum((Zref*BgRef).*Zref,2);
    if any(~isfinite(QrefCell{g})) || any(QrefCell{g} < 0)
        asympt.reason = 'A projected adaptive reference radius is nonfinite.';
        return
    end

    Bembed = zeros(p,p);
    Bembed(obs,obs) = BgRef;
    bgMat(:,g) = Dp'*Bembed(:);

    % TEM-side conditional reconstruction matrices. The empirical sandwich
    % uses the fitted TEM geometry because it is the root of the sample TEM
    % estimating equations.
    SgTEM = STEM(obs,obs);
    SgTEM = (SgTEM+SgTEM')/2;
    if rcond(SgTEM) <= 1e-12
        asympt.reason = ['A fitted TEM pattern scatter block is numerically ' ...
            'singular in the adaptive calibration.'];
        return
    end
    BgTEM = SgTEM\eye(pg);
    BgTEM = (BgTEM+BgTEM')/2;

    Lg = zeros(p,pg);
    Lg(obs,:) = eye(pg);
    Vg = zeros(p,p);
    if ~isempty(mis)
        Ag = STEM(mis,obs)/SgTEM;
        Cg = STEM(mis,mis)-Ag*STEM(obs,mis);
        Cg = (Cg+Cg')/2;
        Lg(mis,:) = Ag;
        Vg(mis,mis) = Cg;
    end

    Lcell{g} = Lg;
    Vcell{g} = Vg;
    Mcell{g} = Lg*SgTEM*Lg';
    BgTEMcell{g} = BgTEM;
    active(g) = true;
end

if ~any(active)
    asympt.reason = 'No active missingness pattern is available at the adaptive threshold.';
    return
end

% ---------------------------------------------------------------------
% Empirical radial factors and their fixed-cutoff influence functions.
% ---------------------------------------------------------------------
if adaptivepool
    dims = unique(pobsDiag(active));
    nGroups = numel(dims);

    dimFD = zeros(nGroups,1);
    nPatternsFD = zeros(nGroups,1);
    aFD = NaN(nGroups,1);
    HFD = NaN(nGroups,1);
    kFD = NaN(nGroups,1);
    fFD = NaN(nGroups,1);
    bwFD = NaN(nGroups,1);
    rhoFD = NaN(nGroups,1);
    nrefFD = zeros(nGroups,1);
    nrefkeptFD = zeros(nGroups,1);

    for jd = 1:nGroups
        kd = dims(jd);
        members = find(active & pobsDiag==kd);
        mk = numel(members);
        ak = local_tem_inv_adjust(cthr,kd,p,n,method);

        % Exactly as in mdTEM with adaptivepool=true, concatenate one
        % projected reference-radius vector for each retained pattern of the
        % same observed dimension. The asymptotic influence below, however,
        % is aggregated by complete row and therefore preserves dependence.
        qpool = vertcat(QrefCell{members});
        insidePool = isfinite(qpool) & qpool <= ak;
        ninside = sum(insidePool);
        minNeeded = max(adaptiveminref,kd+1);
        if ninside < minNeeded
            asympt.reason = sprintf(['Pooled adaptive analytical calibration ' ...
                'has only %d projected reference radii inside the cutoff for ' ...
                'observed dimension %d; at least adaptiveminref=%d are required.'], ...
                ninside,kd,minNeeded);
            return
        end

        Hk = ninside/numel(qpool);
        kk = mean(qpool(insidePool))/kd;
        [fk,bwk] = local_radial_density_reflect(qpool,ak);
        if ~isfinite(Hk) || Hk <= 0 || ~isfinite(kk) || kk <= 0 || ...
                ~isfinite(fk) || fk <= 0
            asympt.reason = ['Unable to estimate a pooled adaptive radial ' ...
                'factor or boundary density for one observed dimension.'];
            return
        end
        rhok = -kk + (ak*fk/Hk)*(ak/kd-kk);

        % Observation-level direct factor influence. Each complete row
        % contributes the sum of its correlated projections; there is no
        % artificial factor mk in the effective sample size.
        directCC = zeros(ncc,1);
        for jj = 1:mk
            g = members(jj);
            qg = QrefCell{g};
            ing = qg <= ak;
            directCC = directCC + ing.*(qg-kd*kk);
        end
        directCC = directCC/(qhat*kd*mk*Hk);

        psiK = zeros(n,1);
        psiK(completeIdx) = directCC;

        % Plug-in robust-scatter contribution uses the average pattern
        % direction because the pooled factor gives equal weight to the mk
        % projected reference samples.
        bbar = mean(bgMat(:,members),2);
        psiK = psiK + (rhok/kd)*(PsiCfull*bbar);

        % One pooled factor perturbs all patterns of this dimension at once.
        factorCoeff = zeros(s,1);
        for jj = 1:mk
            g = members(jj);
            factorCoeff = factorCoeff + ...
                piDiag(g)*(Hk/kk)*local_vech(Mcell{g});

            HDiag(g) = Hk;
            kDiag(g) = kk;
            fDiag(g) = fk;
            bwDiag(g) = bwk;
            rhoDiag(g) = rhok;
            nrefDiag(g) = numel(qpool);
            nrefkeptDiag(g) = ninside;
            poolSizeDiag(g) = mk;
            factorGroupDiag(g) = jd;
            kVec(g) = kk;
        end
        XiFactor = XiFactor - psiK*factorCoeff';

        dimFD(jd) = kd;
        nPatternsFD(jd) = mk;
        aFD(jd) = ak;
        HFD(jd) = Hk;
        kFD(jd) = kk;
        fFD(jd) = fk;
        bwFD(jd) = bwk;
        rhoFD(jd) = rhok;
        nrefFD(jd) = numel(qpool);
        nrefkeptFD(jd) = ninside;
    end

    asympt.TEM.nFactorGroups = nGroups;
    asympt.TEM.factorDiagnostics = table(dimFD,nPatternsFD,aFD,HFD,kFD, ...
        fFD,bwFD,rhoFD,nrefFD,nrefkeptFD,'VariableNames', ...
        {'pobs','nPatterns','athr','H','kappa','density','bandwidth', ...
        'rho','nref','nrefkept'});
else
    members = find(active);
    nGroups = numel(members);
    for jj = 1:nGroups
        g = members(jj);
        pg = pobsDiag(g);
        ag = aDiag(g);
        Qref = QrefCell{g};
        inside = Qref <= ag;
        ninside = sum(inside);
        minNeeded = max(adaptiveminref,pg+1);
        if ninside < minNeeded
            asympt.reason = sprintf(['Adaptive analytical calibration has only %d ' ...
                'reference radii inside the cutoff for an active p_g=%d pattern; ' ...
                'at least adaptiveminref=%d are required.'], ...
                ninside,pg,minNeeded);
            return
        end

        Hg = ninside/ncc;
        kg = mean(Qref(inside))/pg;
        [fg,bwg] = local_radial_density_reflect(Qref,ag);
        if ~isfinite(Hg) || Hg <= 0 || ~isfinite(kg) || kg <= 0 || ...
                ~isfinite(fg) || fg <= 0
            asympt.reason = ['Unable to estimate an adaptive radial factor or ' ...
                'boundary density for one active pattern.'];
            return
        end
        rho = -kg + (ag*fg/Hg)*(ag/pg-kg);

        psiKg = zeros(n,1);
        ccRad = zeros(ncc,1);
        ccRad(inside) = (Qref(inside)-pg*kg)/(qhat*pg*Hg);
        psiKg(completeIdx) = ccRad;
        psiKg = psiKg + (rho/pg)*(PsiCfull*bgMat(:,g));

        factorCoeff = piDiag(g)*(Hg/kg)*local_vech(Mcell{g});
        XiFactor = XiFactor - psiKg*factorCoeff';

        HDiag(g) = Hg;
        kDiag(g) = kg;
        fDiag(g) = fg;
        bwDiag(g) = bwg;
        rhoDiag(g) = rho;
        nrefDiag(g) = ncc;
        nrefkeptDiag(g) = ninside;
        poolSizeDiag(g) = 1;
        factorGroupDiag(g) = jj;
        kVec(g) = kg;
    end
    asympt.TEM.nFactorGroups = nGroups;
end

% Effective adaptive Jacobian. The structural form is unchanged by pooling:
% patterns are still summed separately through pi_g, L_g, M_g and B_g,
% while H, kappa and the boundary density are common within a pooled group.
for g = find(active)'
    pg = pobsDiag(g);
    ag = aDiag(g);
    Hg = HDiag(g);
    kg = kDiag(g);
    fg = fDiag(g);
    obs = obsCell{g};
    Lg = Lcell{g};
    Mg = Mcell{g};
    BgTEM = BgTEMcell{g};

    % Location Jacobian under the same local shell conditions used by the
    % adaptive TEM theorem.
    Pg = zeros(pg,p);
    Pg(:,obs) = eye(pg);
    jmuCoeff = -Hg + 2*ag*fg/pg;
    Jmu = Jmu + piDiag(g)*jmuCoeff*(Lg*Pg);

    ug = -Hg + 2*fg*ag^2/(kg*pg*(pg+2));
    vg = fg*(ag^2/(kg*pg*(pg+2))-ag/pg);
    for h = 1:s
        Hmat = reshape(Dp(:,h),p,p);
        Hsub = Hmat(obs,obs);
        tg = trace(BgTEM*Hsub);
        JH = ug*(Lg*Hsub*Lg') + vg*tg*Mg;
        J(:,h) = J(:,h)+piDiag(g)*local_vech(JH);
    end
end

asympt.TEM.nActivePatterns = sum(active);
asympt.TEM.patternDiagnostics = table(pobsDiag,piDiag,nkeptDiag,aDiag,HDiag, ...
    kDiag,fDiag,bwDiag,rhoDiag,nrefDiag,nrefkeptDiag,poolSizeDiag, ...
    factorGroupDiag,'VariableNames', ...
    {'pobs','pi','nkept','athr','H','kappa','density','bandwidth','rho', ...
    'nref','nrefkept','poolSize','factorGroup'});

Jrcond = rcond(J);
asympt.TEM.Jrcond = Jrcond;
if ~isfinite(Jrcond) || Jrcond <= 1e-12
    asympt.reason = ['The adaptive effective TEM scatter Jacobian is ' ...
        'numerically singular; analytical calibration not computed.'];
    return
end
Jmurcond = rcond(Jmu);
asympt.TEM.Jmurcond = Jmurcond;
jointLocationAvailable = isfinite(Jmurcond) && Jmurcond > 1e-12;
if ~jointLocationAvailable
    asympt.jointReason = ['The adaptive effective TEM location Jacobian is ' ...
        'numerically singular; generic joint MDP-I calibration is unavailable, ' ...
        'while the scatter-only analytical results remain usable.'];
end

% Base innovation evaluated at the fitted TEM root. At the population
% common target this is asymptotically equivalent to using the complete-case
% plug-in, but the TEM plug-in respects the sample estimating equation and
% should have an empirical mean close to zero before explicit centering.
for i = 1:n
    if ~w(i)
        continue
    end
    g = ic(i);
    if ~active(g)
        continue
    end
    obs = obsCell{g};
    zi = Y(i,obs)'-muTEM(obs);
    xhat = Lcell{g}*zi;
    XiMu(i,:) = xhat';
    Ui = (xhat*xhat')/kVec(g)+Vcell{g};
    XiBase(i,:) = local_vech(Ui-STEM)';
end
XiEff = XiBase+XiFactor;
asympt.TEM.estimatingEquationMeanNorm = norm(mean(XiEff,1));

% Joint location--scatter estimator-discrepancy influence. The scatter block
% remains the covariance used by the Hausman omnibus statistic.
PsiTEMfull = -(J\XiEff')';
PsiDelta = PsiTEMfull-PsiCfull;
meanPsiDelta = mean(PsiDelta,1);
PsiDeltaCentered = PsiDelta-meanPsiDelta;
OmegaDelta = (PsiDeltaCentered'*PsiDeltaCentered)/n;
OmegaDelta = (OmegaDelta+OmegaDelta')/2;
asympt.OmegaDelta = OmegaDelta;
asympt.discrepancyIFMeanNorm = norm(meanPsiDelta);
if jointLocationAvailable
    PsiMuTEMfull = -(Jmu\XiMu')';
    PsiDeltaMu = PsiMuTEMfull-PsiMuCfull;
    meanPsiDeltaMu = mean(PsiDeltaMu,1);
    PsiDeltaMuCentered = PsiDeltaMu-meanPsiDeltaMu;
    OmegaDeltaMu = (PsiDeltaMuCentered'*PsiDeltaMuCentered)/n;
    OmegaDeltaMuSigma = (PsiDeltaMuCentered'*PsiDeltaCentered)/n;
    OmegaJoint = [OmegaDeltaMu OmegaDeltaMuSigma; ...
        OmegaDeltaMuSigma' OmegaDelta];
    OmegaJoint = (OmegaJoint+OmegaJoint')/2;
    asympt.OmegaDeltaMu = OmegaDeltaMu;
    asympt.OmegaDeltaMuSigma = OmegaDeltaMuSigma;
    asympt.OmegaJoint = OmegaJoint;
    asympt.locationDiscrepancyIFMeanNorm = norm(meanPsiDeltaMu);
end

% Project the adaptive TEM scatter influence onto the ordinary MDP-U
% coefficient. Since psi_TEM=-J_eff^{-1}xi_eff and TD=-a0'(psi_TEM-psi_C),
% the TEM contribution is b'xi_eff, where J_eff' b=a0.
b = J'\a0;
temBase = XiBase*b;
temFactor = XiFactor*b;
temContribution = temBase+temFactor;
meanBaseRaw = mean(temBase);
meanFactorRaw = mean(temFactor);
meanTemRaw = mean(temContribution);
temBase = temBase-meanBaseRaw;
temFactor = temFactor-meanFactorRaw;
temContribution = temContribution-meanTemRaw;

% Robust complete-case contribution to TD,mean.
ccContribution = PsiCfull*a0;
meanCCRaw = mean(ccContribution);
ccContribution = ccContribution-meanCCRaw;

zetaD = temContribution+ccContribution;
meanTotalRaw = meanTemRaw+meanCCRaw;
zetaD = zetaD-mean(zetaD);

sigmaBase2 = mean(temBase.^2);
sigmaFactor2 = mean(temFactor.^2);
crossBaseFactor = mean(temBase.*temFactor);
sigmaTEM2 = mean(temContribution.^2);
sigmaCC2 = mean(ccContribution.^2);
crossTEMCC = mean(temContribution.*ccContribution);
sigmaD2 = mean(zetaD.^2);
sigmaD2FromOmega = a0'*OmegaDelta*a0;
asympt.omnibusVarianceIdentityResidual = sigmaD2FromOmega-sigmaD2;

asympt.TEM.available = true;
asympt.TEM.sigma2 = sigmaTEM2;
asympt.TEM.sigma2Base = sigmaBase2;
asympt.TEM.sigma2Factor = sigmaFactor2;
asympt.TEM.twiceBaseFactor = 2*crossBaseFactor;
asympt.TEM.meanInfluenceBeforeCentering = meanTemRaw;
asympt.TEM.meanBaseBeforeCentering = meanBaseRaw;
asympt.TEM.meanFactorBeforeCentering = meanFactorRaw;
asympt.completeCase.sigma2 = sigmaCC2;
asympt.completeCase.crossTEMCC = crossTEMCC;
asympt.completeCase.meanInfluenceBeforeCentering = meanCCRaw;
asympt.meanInfluenceBeforeCentering = meanTotalRaw;

if adaptivepool
    poolText = [' Dimension-pooled radial factors are used; their influence ' ...
        'is evaluated at the complete-observation level so dependence among ' ...
        'projections of the same complete row is retained.'];
else
    poolText = ' Pattern-specific adaptive radial factors are used.';
end

methodText = [' Parameter-free distance mapping method=''' method ''' is ' ...
    'handled through its inverse raw-distance cutoff.'];
if strcmpi(robustClass,'FS')
    asympt.theoryStatus = ['Adaptive TEM sandwich under the local shell-moment conditions.' ...
        methodText ...
        poolText ' The complete-case FS scatter influence is evaluated ' ...
        'analytically conditional on the full-sample final FS classification; ' ...
        'the data-dependent FS stopping/classification rule itself is not ' ...
        'differentiated and remains outside the fixed-fraction FS theorem.'];
else
    asympt.theoryStatus = ['Adaptive TEM sandwich under the local shell-moment conditions.' ...
        methodText poolText ' The empirical radial-factor influence and robust ' ...
        'complete-case scatter influence are both included.'];
end

if ~isfinite(sigmaD2) || sigmaD2 < 0
    asympt.reason = 'The adaptive end-to-end sandwich variance is not finite.';
    return
end

tolVar = 1e-10*max(1,sigmaTEM2+sigmaCC2);
asympt.available = true;
asympt.reason = '';
asympt.qhat = qhat;
asympt.degenerate = sigmaD2 <= tolVar;
asympt.baseSigmaD2 = sigmaD2;
asympt.baseKappa = baseKappa;
asympt.radialMeanQ = meanQ0;
asympt.invTauHat = invTauHat;

asympt.TDmean.sigma2 = sigmaD2;
asympt.TDmean.se = sqrt(sigmaD2/n);
asympt.TDmean.components = struct('TEM',sigmaTEM2, ...
    'completeCase',sigmaCC2,'twiceCross',2*crossTEMCC, ...
    'TEMbase',sigmaBase2,'TEMfactor',sigmaFactor2, ...
    'TEMtwiceBaseFactor',2*crossBaseFactor);

asympt.TLmean.kappa = baseKappa;
asympt.TLmean.sigma2 = baseKappa^2*sigmaD2;
asympt.TLmean.se = sqrt(asympt.TLmean.sigma2/n);

if asympt.degenerate
    asympt.TDmean.z = NaN;
    asympt.TDmean.pvalue = NaN;
    asympt.TLmean.z = NaN;
    asympt.TLmean.pvalue = NaN;
else
    asympt.TDmean.z = Tobs(4)/asympt.TDmean.se;
    asympt.TDmean.pvalue = erfc(abs(asympt.TDmean.z)/sqrt(2));
    asympt.TLmean.z = Tobs(2)/asympt.TLmean.se;
    asympt.TLmean.pvalue = erfc(abs(asympt.TLmean.z)/sqrt(2));
end
end

% -------------------------------------------------------------------------
function [f,h] = local_radial_density_reflect(q,a)
%local_radial_density_reflect Boundary density estimate on [0,Inf).
%
% A Gaussian kernel with reflection at zero is used:
%   f(a)=mean{phi((a-Q)/h)+phi((a+Q)/h)}/h.
% The bandwidth is Silverman's robust rule with deterministic fallbacks.

q = q(:);
q = q(isfinite(q) & q >= 0);
f = NaN;
h = NaN;
if isempty(q) || ~isfinite(a) || a < 0
    return
end

nq = numel(q);
sq = std(q,1);
if nq >= 4
    qq = quantile(q,[0.25 0.75]);
    riqr = (qq(2)-qq(1))/1.349;
else
    riqr = NaN;
end
scale = sq;
if isfinite(riqr) && riqr > 0
    if isfinite(scale) && scale > 0
        scale = min(scale,riqr);
    else
        scale = riqr;
    end
end
if ~isfinite(scale) || scale <= 0
    scale = max(median(abs(q-median(q))),sqrt(eps)*max(1,a));
end
h = 0.9*scale*nq^(-1/5);
h = max(h,sqrt(eps)*max(1,a));

u1 = (a-q)/h;
u2 = (a+q)/h;
phi1 = exp(-0.5*u1.^2)/sqrt(2*pi);
phi2 = exp(-0.5*u2.^2)/sqrt(2*pi);
f = mean(phi1+phi2)/h;
end

% -------------------------------------------------------------------------
function [PsiMuC,PsiC,info,reason] = local_cc_joint_influence(Ycc,alpha, ...
    robustClass,robustEff,robustBonflev,outRob)
%local_cc_joint_influence Robust complete-case location/scatter influence.
%
% For FS, use an analytical frozen-classification perturbation of the final
% retained mean and covariance. For MCD and MM, a single delete-one loop is
% used to obtain both location and scatter influence values.

info = struct('influenceMethod','', ...
    'nReferenceInliers',NaN,'nReferenceOutliers',NaN, ...
    'frozenClassification',false,'relCovFrozenVsFS',NaN, ...
    'relLocFrozenVsFS',NaN);
PsiMuC = [];
PsiC = [];
reason = '';

if strcmpi(robustClass,'FS')
    [PsiMuC,PsiC,info,reason] = local_cc_joint_fs_frozen(Ycc,outRob);
else
    [PsiMuC,PsiC,reason] = local_cc_joint_jackknife(Ycc,alpha,robustClass, ...
        robustEff,robustBonflev);
    if isempty(reason)
        info.influenceMethod = 'delete-one jackknife';
        info.frozenClassification = false;
    end
end
end

% -------------------------------------------------------------------------
function [PsiMuC,PsiC,info,reason] = local_cc_joint_fs_frozen(Ycc,outRob)
%local_cc_joint_fs_frozen Analytical frozen-classification FS influence.
%
% Conditional on the final inlier set I, h=|I|, the location influence is
%
%   Psi_mu,i = (m-1)/(h-1) * (y_i-ybar_I),  i in I,
%              0,                            i not in I.
%
% The scatter influence is the corresponding exact frozen-classification
% delete-one perturbation used previously by mdMDPtest.

m = size(Ycc,1);
p = size(Ycc,2);
s = p*(p+1)/2;
PsiMuC = [];
PsiC = [];
reason = '';
info = struct('influenceMethod','FS analytic frozen-classification', ...
    'nReferenceInliers',NaN,'nReferenceOutliers',NaN, ...
    'frozenClassification',true,'relCovFrozenVsFS',NaN, ...
    'relLocFrozenVsFS',NaN);

if m < 5
    reason = 'Too few complete observations for the FS joint influence.';
    return
end
if isempty(outRob) || ~isstruct(outRob) || ~isfield(outRob,'cov') || ...
        ~isfield(outRob,'loc')
    reason = ['The full-sample FS fit is required for the analytical ' ...
        'frozen-classification joint influence.'];
    return
end

outMask = local_cc_outlier_mask(outRob,m);
inMask = ~outMask;
h = sum(inMask);
info.nReferenceInliers = h;
info.nReferenceOutliers = sum(outMask);
if h < 3
    reason = ['Too few observations remain in the final FS inlier set for ' ...
        'the frozen-classification influence.'];
    return
end

Yin = Ycc(inMask,:);
muI = mean(Yin,1)';
SI = cov(Yin);
SI = (SI+SI')/2;
muFS = outRob.loc(:);
SFS = (outRob.cov+outRob.cov')/2;
if any(~isfinite(SI(:))) || any(~isfinite(SFS(:))) || ...
        any(~isfinite(muI)) || any(~isfinite(muFS))
    reason = 'The final FS location or scatter contains nonfinite values.';
    return
end

info.relLocFrozenVsFS = norm(muI-muFS)/max(norm(muFS),1);
info.relCovFrozenVsFS = norm(SI-SFS,'fro')/max(norm(SFS,'fro'),eps);
if info.relLocFrozenVsFS > 1e-8 || info.relCovFrozenVsFS > 1e-8
    reason = ['FSM location/scatter differs from the mean/covariance of the ' ...
        'final FS inlier set; the analytical frozen-classification influence ' ...
        'is therefore not applied.'];
    return
end

PsiMuC = zeros(m,p);
PsiC = zeros(m,s);
vS = local_vech(SI);
inIdx = find(inMask);
coefMu = (m-1)/(h-1);
coefS = (m-1)/(h-2);
radialCoef = h/(h-1);
for j = 1:numel(inIdx)
    i = inIdx(j);
    x = Ycc(i,:)'-muI;
    PsiMuC(i,:) = (coefMu*x)';
    vA = local_vech(x*x');
    PsiC(i,:) = (coefS*(radialCoef*vA-vS))';
end
PsiMuC = PsiMuC-mean(PsiMuC,1);
PsiC = PsiC-mean(PsiC,1);
if any(~isfinite(PsiMuC(:))) || any(~isfinite(PsiC(:)))
    reason = 'The analytical frozen-FS joint influence is nonfinite.';
    PsiMuC = [];
    PsiC = [];
end
end

% -------------------------------------------------------------------------
function [PsiMuC,PsiC,reason] = local_cc_joint_jackknife(Ycc,alpha,robustClass, ...
    robustEff,robustBonflev)
%local_cc_joint_jackknife Delete-one robust location/scatter influence.
%
% The same delete-one fits supply both location and scatter perturbations,
% avoiding two separate robust jackknife loops.

ncc = size(Ycc,1);
p = size(Ycc,2);
s = p*(p+1)/2;
PsiMuC = [];
PsiC = [];
reason = '';
if ncc < 5
    reason = 'Too few complete observations for the robust delete-one influence.';
    return
end

looMu = NaN(ncc,p);
looS = NaN(ncc,s);
rngState = rng;
cleanupObj = onCleanup(@() rng(rngState)); %#ok<NASGU>
for i = 1:ncc
    keep = true(ncc,1);
    keep(i) = false;
    try
        rng(rngState)
        [muMinus,SigMinus] = local_complete_case_fit(Ycc(keep,:),alpha, ...
            robustClass,robustEff,robustBonflev);
        muMinus = muMinus(:);
        if any(~isfinite(muMinus)) || any(~isfinite(SigMinus(:)))
            reason = sprintf(['Nonfinite robust location/scatter in delete-one ' ...
                'fit for complete observation %d.'],i);
            return
        end
        looMu(i,:) = muMinus';
        looS(i,:) = local_vech((SigMinus+SigMinus')/2)';
    catch ME
        reason = sprintf(['Robust delete-one fit failed for complete ' ...
            'observation %d: %s'],i,ME.message);
        return
    end
end
if any(~isfinite(looMu(:))) || any(~isfinite(looS(:)))
    reason = 'The robust delete-one joint influence contains nonfinite values.';
    return
end
PsiMuC = -(ncc-1)*(looMu-mean(looMu,1));
PsiC = -(ncc-1)*(looS-mean(looS,1));
end

% -------------------------------------------------------------------------
function asympt = local_tem_asymptotic(Y, maskMiss, completeIdx, outFit, ...
    alpha, method, robustClass, robustEff, robustBonflev, Tobs, eps0, ...
    outRob, frozenEligibility)
%local_tem_asymptotic Analytical alpha>0 TEM sandwich benchmark.
%
% The implementation follows the scalar influence-function representation
%
%   J_eff' * b = aSigma,
%   -aSigma' * psi_TEM(W) = b' * xi(W),
%
% where xi(W)=vech{w_R(U_R-Sigma)}. The complete-case contribution to the
% MDP statistic is C*zeta_C/q. For FS the complete-case scatter influence is
% evaluated analytically conditional on the final full-sample FS
% classification. For MCD and MM the current delete-one jackknife is
% retained.
%
% The Jacobian implementation covers the parameter-free adjusted-distance
% mappings 'pri', 'expScale', 'zMap', 'chiMap' and 'betaMap'. The selected
% mapping enters only through the inverse raw-distance cutoff a_g(c). At the
% correctly consistency-adjusted population target the estimating equation is
% zero for every admissible cutoff, so the first-order contribution from the
% random common cutoff cancels with the cutoff-induced change of the Tallis
% factor. The scatter-dependent 'detMap' mapping would require additional
% derivatives and is deliberately excluded.

[n,p] = size(Y);
qhat = sum(completeIdx)/n;
s = p*(p+1)/2;

if nargin < 13 || isempty(frozenEligibility)
    frozenEligibility = true(n,1);
else
    frozenEligibility = logical(frozenEligibility(:));
    if numel(frozenEligibility) ~= n
        error('FSDA:mdMDPtest:WrongCoupledMask', ...
            'The frozen eligibility mask must have one entry for every row of Y.');
    end
end
frozenMode = any(~frozenEligibility);

asympt = local_unavailable_asymptotic(qhat);
asympt.reason = '';
asympt.VA = NaN;
asympt.mode = 'TEM influence-function sandwich';
asympt.theoryStatus = 'available';
asympt.TEM = struct('available',false,'sigma2',NaN,'threshold',NaN, ...
    'Jrcond',NaN,'Jmurcond',NaN,'nActivePatterns',0,'nPatterns',0);
asympt.completeCase = struct('influenceMethod','', ...
    'sigma2',NaN,'crossTEMCC',NaN,'nComplete',sum(completeIdx), ...
    'nReferenceInliers',NaN,'nReferenceOutliers',NaN, ...
    'frozenClassification',false,'robustBonflev',robustBonflev, ...
    'relCovFrozenVsFS',NaN,'relLocFrozenVsFS',NaN);

supportedMethods = {'pri','expScale','zMap','chiMap','betaMap'};
if strcmpi(method,'detMap')
    asympt.reason = ['For alpha>0 the analytical TEM sandwich deliberately ' ...
        'excludes method=''detMap'' because its adjusted-distance map ' ...
        'depends explicitly on the scatter matrix.'];
    return
elseif ~any(strcmpi(method,supportedMethods))
    asympt.reason = ['For alpha>0 the analytical TEM sandwich is implemented ' ...
        'for ''pri'', ''expScale'', ''zMap'', ''chiMap'' and ''betaMap''.'];
    return
end

if isempty(outFit) || ~isfield(outFit,'weights') || ...
        ~isfield(outFit,'loc') || ~isfield(outFit,'cov')
    asympt.reason = 'The fitted TEM object does not contain the required fields.';
    return
end

muHat = outFit.loc(:);
SigmaHat = (outFit.cov + outFit.cov')/2;
if any(~isfinite(SigmaHat(:))) || rcond(SigmaHat) <= 1e-12
    asympt.reason = 'The fitted TEM scatter matrix is numerically singular.';
    return
end

% Duplication matrix and derivative of log|Sigma|.
Dp = local_duplication_matrix(p);
K = SigmaHat\eye(p);
K = (K+K')/2;
aSigma = Dp' * K(:);

% Reconstruct the common adjusted-distance threshold used by TEM. Current
% mdTEM objects store cthr and the final adjusted distances explicitly. The
% adjusted-distance fallback also makes the analytical code independent of
% the particular parameter-free mapping.
cthr = NaN;
if isfield(outFit,'cthr') && isscalar(outFit.cthr) && isfinite(outFit.cthr)
    cthr = outFit.cthr;
end

if isfield(outFit,'adjustedD2') && numel(outFit.adjustedD2)==n
    d2adj = outFit.adjustedD2(:);
else
    d2adj = local_adjusted_distances(Y,muHat,SigmaHat,method);
end
w = outFit.weights(:) > 0;

if ~isfinite(cthr)
    if ~any(w & isfinite(d2adj))
        asympt.reason = 'Unable to reconstruct the final TEM trimming threshold.';
        return
    end
    cthr = max(d2adj(w & isfinite(d2adj)));
end

% Missingness patterns and empirical pattern probabilities.
[patt,~,ic] = unique(maskMiss,'rows');
G = size(patt,1);
asympt.TEM.nPatterns = G;

J = zeros(s,s);
Jmu = zeros(p,p);
active = false(G,1);
kgVec = NaN(G,1);
Lcell = cell(G,1);
Vcell = cell(G,1);
obsCell = cell(G,1);

% Build the effective scatter Jacobian pattern by pattern. For every
% supported parameter-free mapping, first invert the common adjusted cutoff c
% to the corresponding raw p_g-dimensional squared-distance cutoff a_g(c).
for g = 1:G
    obs = find(~patt(g,:));
    mis = find(patt(g,:));
    pg = numel(obs);
    pig = sum((ic==g) & frozenEligibility)/n;
    obsCell{g} = obs;

    if pg == 0 || pig <= 0
        continue
    end

    ag = local_tem_inv_adjust(cthr,pg,p,n,method);
    if ~isfinite(ag) || ag <= 0
        continue
    end

    gammag = chi2cdf(ag,pg);
    fg = chi2pdf(ag,pg);
    kg = chi2cdf(ag,pg+2)/gammag;

    if ~isfinite(kg) || kg <= 0 || ~isfinite(fg)
        continue
    end

    SigmaG = SigmaHat(obs,obs);
    SigmaG = (SigmaG+SigmaG')/2;
    if rcond(SigmaG) <= 1e-12
        asympt.reason = ['A pattern-specific covariance block is numerically ' ...
            'singular; TEM analytical benchmark not computed.'];
        return
    end
    BG = SigmaG\eye(pg);
    BG = (BG+BG')/2;

    Lg = zeros(p,pg);
    Lg(obs,:) = eye(pg);
    Vg = zeros(p,p);

    if ~isempty(mis)
        Ag = SigmaHat(mis,obs)/SigmaG;
        Cg = SigmaHat(mis,mis)-Ag*SigmaHat(obs,mis);
        Cg = (Cg+Cg')/2;
        Lg(mis,:) = Ag;
        Vg(mis,mis) = Cg;
    end

    % Location Jacobian: J_mumu u = sum_g pi_g
    % (-gamma_g+2*a_g*f_g/p_g) L_g P_g u.
    Pg = zeros(pg,p);
    Pg(:,obs) = eye(pg);
    jmuCoeff = -gammag + 2*ag*fg/pg;
    Jmu = Jmu + pig*jmuCoeff*(Lg*Pg);

    ug = -gammag + 2*ag^2*fg/(pg*(pg+2)*kg);
    vg = fg*(ag^2/(pg*(pg+2)*kg)-ag/pg);
    LSigmaL = Lg*SigmaG*Lg';

    for h = 1:s
        H = reshape(Dp(:,h),p,p);
        Hg = H(obs,obs);
        tg = trace(BG*Hg);
        JH = ug*(Lg*Hg*Lg') + vg*tg*LSigmaL;
        J(:,h) = J(:,h) + pig*local_vech(JH);
    end

    active(g) = true;
    kgVec(g) = kg;
    Lcell{g} = Lg;
    Vcell{g} = Vg;
end

asympt.TEM.threshold = cthr;
asympt.TEM.nActivePatterns = sum(active);

if ~any(active)
    asympt.reason = 'No active missingness pattern is available at the TEM threshold.';
    return
end

Jrcond = rcond(J);
asympt.TEM.Jrcond = Jrcond;
if ~isfinite(Jrcond) || Jrcond <= 1e-12
    asympt.reason = ['The effective TEM scatter Jacobian is numerically ' ...
        'singular; analytical benchmark not computed.'];
    return
end
Jmurcond = rcond(Jmu);
asympt.TEM.Jmurcond = Jmurcond;
jointLocationAvailable = isfinite(Jmurcond) && Jmurcond > 1e-12;
if ~jointLocationAvailable
    asympt.jointReason = ['The effective TEM location Jacobian is numerically ' ...
        'singular; generic joint MDP-I calibration is unavailable, while the ' ...
        'scatter-only analytical results remain usable.'];
end

% b gives the scalar projection needed by TD,mean without constructing the
% full TEM covariance matrix.
b = J'\aSigma;

% Empirical values of xi_i=vech{w_i(U_i-Sigma)}. The theoretical exact
% pattern factor is used here, rather than the finite-sample shrinkage of
% unstable factors that mdTEM may apply internally.
Xi = zeros(n,s);
XiMu = zeros(n,p);
for i = 1:n
    if ~w(i)
        continue
    end
    g = ic(i);
    if ~active(g)
        continue
    end

    obs = obsCell{g};
    zi = Y(i,obs)' - muHat(obs);
    xhat = Lcell{g}*zi;
    XiMu(i,:) = xhat';
    Ui = (xhat*xhat')/kgVec(g) + Vcell{g};
    Xi(i,:) = local_vech(Ui-SigmaHat)';
end

% b'xi is the TEM contribution to the influence function of TD,mean,
% because TD,mean has the opposite sign of aSigma'psi_TEM.
temContribution = Xi*b;
temContribution = temContribution-mean(temContribution);
sigmaTEM2 = mean(temContribution.^2);
asympt.TEM.available = true;
asympt.TEM.sigma2 = sigmaTEM2;

% Complete-case scatter influence. For FS use the analytical frozen-
% classification perturbation; MCD and MM retain the delete-one jackknife.
Ycc = Y(completeIdx,:);
[PsiMuCcc,PsiCcc,ccInfo,ccReason] = local_cc_joint_influence(Ycc,alpha, ...
    robustClass,robustEff,robustBonflev,outRob);
if ~isempty(ccReason)
    asympt.reason = ccReason;
    return
end
zetaC = PsiCcc*aSigma;
asympt.completeCase.influenceMethod = ccInfo.influenceMethod;
asympt.completeCase.nReferenceInliers = ccInfo.nReferenceInliers;
asympt.completeCase.nReferenceOutliers = ccInfo.nReferenceOutliers;
asympt.completeCase.frozenClassification = ccInfo.frozenClassification;
asympt.completeCase.relCovFrozenVsFS = ccInfo.relCovFrozenVsFS;
asympt.completeCase.relLocFrozenVsFS = ccInfo.relLocFrozenVsFS;

PsiCfull = zeros(n,s);
PsiCfull(completeIdx,:) = PsiCcc/qhat;
PsiMuCfull = zeros(n,p);
PsiMuCfull(completeIdx,:) = PsiMuCcc/qhat;

% Joint location--scatter estimator-discrepancy influence. The scatter block
% is also the covariance used by the Hausman omnibus statistic.
PsiTEMfull = -(J\Xi')';
PsiDelta = PsiTEMfull-PsiCfull;
meanPsiDelta = mean(PsiDelta,1);
PsiDeltaCentered = PsiDelta-meanPsiDelta;
OmegaDelta = (PsiDeltaCentered'*PsiDeltaCentered)/n;
OmegaDelta = (OmegaDelta+OmegaDelta')/2;
asympt.OmegaDelta = OmegaDelta;
asympt.discrepancyIFMeanNorm = norm(meanPsiDelta);
if jointLocationAvailable
    PsiMuTEMfull = -(Jmu\XiMu')';
    PsiDeltaMu = PsiMuTEMfull-PsiMuCfull;
    meanPsiDeltaMu = mean(PsiDeltaMu,1);
    PsiDeltaMuCentered = PsiDeltaMu-meanPsiDeltaMu;
    OmegaDeltaMu = (PsiDeltaMuCentered'*PsiDeltaMuCentered)/n;
    OmegaDeltaMuSigma = (PsiDeltaMuCentered'*PsiDeltaCentered)/n;
    OmegaJoint = [OmegaDeltaMu OmegaDeltaMuSigma; ...
        OmegaDeltaMuSigma' OmegaDelta];
    OmegaJoint = (OmegaJoint+OmegaJoint')/2;
    asympt.OmegaDeltaMu = OmegaDeltaMu;
    asympt.OmegaDeltaMuSigma = OmegaDeltaMuSigma;
    asympt.OmegaJoint = OmegaJoint;
    asympt.locationDiscrepancyIFMeanNorm = norm(meanPsiDeltaMu);
end

ccContribution = PsiCfull*aSigma;
ccContribution = ccContribution-mean(ccContribution);

% End-to-end influence of TD,mean:
%   zeta_D = b'xi + C*zeta_C/q.
zetaD = temContribution + ccContribution;
zetaD = zetaD-mean(zetaD);

sigmaCC2 = mean(ccContribution.^2);
crossTEMCC = mean(temContribution.*ccContribution);
sigmaD2 = mean(zetaD.^2);
sigmaD2FromOmega = aSigma'*OmegaDelta*aSigma;
asympt.omnibusVarianceIdentityResidual = sigmaD2FromOmega-sigmaD2;

asympt.completeCase.sigma2 = sigmaCC2;
asympt.completeCase.crossTEMCC = crossTEMCC;
asympt.completeCase.nComplete = size(Ycc,1);

methodText = [' The parameter-free distance mapping method=''' method ''' is ' ...
    'handled through its inverse raw-distance cutoff.'];
if strcmpi(robustClass,'FS')
    asympt.theoryStatus = ['TEM influence calculation is analytical.' ...
        methodText ' The ' ...
        'complete-case FS scatter influence is evaluated analytically ' ...
        'conditional on the full-sample final FS classification; the ' ...
        'data-dependent FS stopping/classification rule itself is not ' ...
        'differentiated.'];
else
    asympt.theoryStatus = ['TEM influence calculation is analytical.' ...
        methodText ' The robust complete-case scatter influence is evaluated ' ...
        'by delete-one ' ...
        'jackknife.'];
end
if frozenMode
    asympt.theoryStatus = [asympt.theoryStatus ' The externally imposed ' ...
        'eligibility mask is frozen in this working sandwich: eligible ' ...
        'pattern masses and realized product weights are used, while the ' ...
        'selection mechanism generating the mask is not differentiated.'];
end

if ~isfinite(sigmaD2) || sigmaD2 < 0
    asympt.reason = 'The end-to-end TEM sandwich variance is not finite.';
    return
end

% Stabilized log-ratio coefficient.
kappa = integral(@(q) local_kappa_integrand(q,p,eps0),0,Inf, ...
    'RelTol',1e-10,'AbsTol',1e-12)/p;

tolVar = 1e-10*max(1,sigmaTEM2+sigmaCC2);
asympt.available = true;
asympt.reason = '';
asympt.qhat = qhat;
asympt.degenerate = sigmaD2 <= tolVar;

asympt.TDmean.sigma2 = sigmaD2;
asympt.TDmean.se = sqrt(sigmaD2/n);
asympt.TDmean.components = struct('TEM',sigmaTEM2, ...
    'completeCase',sigmaCC2,'twiceCross',2*crossTEMCC);

asympt.TLmean.kappa = kappa;
asympt.TLmean.sigma2 = kappa^2*sigmaD2;
asympt.TLmean.se = sqrt(asympt.TLmean.sigma2/n);

if asympt.degenerate
    asympt.TDmean.z = NaN;
    asympt.TDmean.pvalue = NaN;
    asympt.TLmean.z = NaN;
    asympt.TLmean.pvalue = NaN;
else
    asympt.TDmean.z = Tobs(4)/asympt.TDmean.se;
    asympt.TDmean.pvalue = erfc(abs(asympt.TDmean.z)/sqrt(2));
    asympt.TLmean.z = Tobs(2)/asympt.TLmean.se;
    asympt.TLmean.pvalue = erfc(abs(asympt.TLmean.z)/sqrt(2));
end
end

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function v = local_vech(S)
%local_vech Lower-triangular half-vectorization, column by column.

p = size(S,1);
v = zeros(p*(p+1)/2,1);
k = 0;
for j = 1:p
    for i = j:p
        k = k+1;
        v(k) = S(i,j);
    end
end
end

% -------------------------------------------------------------------------
function asympt = local_classical_asymptotic(Y, maskMiss, completeIdx, ...
    muHat, SigmaHat, nComplete, Tobs, TDordinary, eps0, d2cc, ...
    preferGeneralF)
%local_classical_asymptotic Analytical benchmark for alpha=0.
%
% For ordinary untrimmed MDP-U, the primary TDmean calibration is the
% feasible general-F sandwich
%
%   sqrt(n)*TDmean -> N(0,sigmaD,F^2),
%
% valid under MCAR and finite fourth moments. The all-available Gaussian
% observed-data score is used as a quasi-score and its empirical sandwich
% covariance estimates the estimator-discrepancy variance. Gaussianity and
% ellipticity are not required for this TDmean result.
%
% The Gaussian-information specialization
%
%   sigmaD,N^2 = 2*p/qhat - aSigma'*(IA\aSigma)
%
% is retained in asympt.gaussian as a diagnostic. Under ellipticity,
% TLmean is a scalar multiple of TDmean to first order. The primary TLmean
% calibration therefore uses the same feasible general-F variance multiplied
% by an empirical radial coefficient estimated from the complete-case
% Mahalanobis distances d2cc. This does not assume Gaussianity.
%
% For aggregation='fixed', preferGeneralF is false. The caller then applies
% the existing Gaussian Tallis fixed-cutoff scaling. The untrimmed general-F
% TDmean result is still returned in asympt.generalF for diagnostics, but is
% not used to calibrate the filtered statistic.

[n,p] = size(Y);
qhat = nComplete/n;
s = p*(p+1)/2;

asympt = local_unavailable_asymptotic(qhat);
asympt.reason = '';
asympt.mode = '';
asympt.primaryCalibration = '';
asympt.radialModel = '';
asympt.theoryStatus = '';
asympt.VA = NaN;

muHat = muHat(:);
SigmaHat = (SigmaHat + SigmaHat')/2;

if numel(muHat) ~= p || ~isequal(size(SigmaHat),[p p]) || ...
        any(~isfinite(muHat)) || any(~isfinite(SigmaHat(:)))
    asympt.reason = 'The alpha=0 all-available fit has invalid dimensions or nonfinite values.';
    return
end

[~,flagSigma] = chol(SigmaHat);
if flagSigma ~= 0
    asympt.reason = 'The alpha=0 all-available covariance matrix is not positive definite.';
    return
end

% Duplication matrix and MDP scatter direction.
Dp = local_duplication_matrix(p);
K = SigmaHat\eye(p);
K = (K + K')/2;
aSigma = Dp' * K(:);

% Missingness patterns, Gaussian quasi-information bread, and observed-data
% scatter quasi-scores evaluated at the common all-available fit.
[patt,~,ic] = unique(maskMiss,'rows','stable');
G = size(patt,1);
patternCounts = accumarray(ic,1,[G 1]);
patternProportions = patternCounts/n;

JA = zeros(s,s);
score = zeros(n,s);
idxFull = reshape(1:p*p,p,p);

for g = 1:G
    rows = find(ic == g);
    ng = numel(rows);
    pig = ng/n;

    obs = ~patt(g,:);
    pg = sum(obs);
    if pg == 0
        % An all-missing row has zero score and zero information.
        continue
    end

    SigmaG = SigmaHat(obs,obs);
    SigmaG = (SigmaG + SigmaG')/2;
    [~,flagG] = chol(SigmaG);
    if flagG ~= 0
        asympt.reason = sprintf(['Observed covariance submatrix is not positive ' ...
            'definite for missingness pattern %d.'],g);
        return
    end

    BG = SigmaG\eye(pg);
    BG = (BG + BG')/2;

    % Dg maps vech(Sigma) to vec(Sigma_g).
    idxG = idxFull(obs,obs);
    Dg = Dp(idxG(:),:);

    Ir = 0.5 * Dg' * kron(BG,BG) * Dg;
    JA = JA + pig*Ir;

    % Pattern-r Gaussian scatter quasi-score:
    % 0.5*Dg'*vec{BG(x_g x_g'-Sigma_g)BG}.
    Xg = Y(rows,obs) - muHat(obs)';
    W = Xg*BG;
    for krow = 1:ng
        M = W(krow,:)'*W(krow,:) - BG;
        score(rows(krow),:) = (0.5*Dg'*M(:))';
    end
end

JA = (JA + JA')/2;
rcondJA = rcond(JA);

% ---------------------------------------------------------------------
% Feasible general-F sandwich for ordinary untrimmed TDmean.
% ---------------------------------------------------------------------
generalF = struct;
generalF.available = false;
generalF.reason = '';
generalF.sigma2 = NaN;
generalF.se = NaN;
generalF.z = NaN;
generalF.pvalue = NaN;
generalF.OmegaDelta = NaN(s,s);
generalF.aSigma = aSigma;
generalF.JA = JA;
generalF.rcondJA = rcondJA;
generalF.scoreMeanNorm = NaN;
generalF.deltaMeanNorm = NaN;
generalF.varianceIdentityResidual = NaN;
generalF.degenerate = false;
generalF.appliesToAggregation = 'ordinary';
generalF.TDmeanObserved = TDordinary;
generalF.theoryStatus = ['Feasible first-order general-F sandwich for ' ...
    'untrimmed ordinary-mean MDP-U under MCAR and finite fourth moments.'];

if ~isfinite(rcondJA) || rcondJA <= 1e-12
    generalF.reason = ['Estimated scatter bread J_A is numerically singular; ' ...
        'general-F calibration not computed.'];
else
    % All-available scatter influence.
    psiA = (JA\score')';

    % Complete-case covariance influence, evaluated at the common
    % all-available quasi-likelihood geometry:
    % C_i/q * vech{(Y_i-mu)(Y_i-mu)'-Sigma}.
    psiC = zeros(n,s);
    rowsC = find(completeIdx);
    XC = Y(rowsC,:) - muHat';

    for krow = 1:numel(rowsC)
        H = XC(krow,:)'*XC(krow,:) - SigmaHat;
        psiC(rowsC(krow),:) = local_vech(H)'/qhat;
    end

    delta = psiA - psiC;
    scoreMean = mean(score,1);
    deltaMean = mean(delta,1);
    deltaCentered = delta - deltaMean;

    OmegaDelta = (deltaCentered'*deltaCentered)/n;
    OmegaDelta = (OmegaDelta + OmegaDelta')/2;

    sigmaD2F = aSigma'*OmegaDelta*aSigma;

    % Independent scalar-IF identity: IF(TDmean)=-aSigma'*delta.
    zeta = -(deltaCentered*aSigma);
    sigmaD2Scalar = mean(zeta.^2);
    varianceIdentityResidual = sigmaD2Scalar - sigmaD2F;

    tolVarF = 1e-10*max(1,trace(OmegaDelta)*max(1,norm(aSigma)^2));

    generalF.OmegaDelta = OmegaDelta;
    generalF.scoreMeanNorm = norm(scoreMean);
    generalF.deltaMeanNorm = norm(deltaMean);
    generalF.varianceIdentityResidual = varianceIdentityResidual;

    if sigmaD2F < 0 && sigmaD2F >= -tolVarF
        sigmaD2F = 0;
    end

    if ~isfinite(sigmaD2F) || sigmaD2F < 0
        generalF.reason = ['Estimated general-F sandwich variance is not ' ...
            'finite and nonnegative.'];
    elseif sigmaD2F <= tolVarF
        generalF.degenerate = true;
        generalF.reason = ['Estimated general-F sandwich variance is ' ...
            'numerically degenerate.'];
    else
        generalF.available = true;
        generalF.reason = '';
        generalF.sigma2 = sigmaD2F;
        generalF.se = sqrt(sigmaD2F/n);
        generalF.z = TDordinary/generalF.se;
        generalF.pvalue = erfc(abs(generalF.z)/sqrt(2));
    end
end

% ---------------------------------------------------------------------
% Elliptical radial extension for ordinary untrimmed TLmean.
% ---------------------------------------------------------------------
ellipticalTL = struct('available',false,'reason','', ...
    'calibration','ellipticalRadialGeneralF','kappa',NaN,'lambda',NaN, ...
    'coefficient',NaN,'sigma2',NaN,'se',NaN,'z',NaN,'pvalue',NaN, ...
    'radialMeanQ',NaN,'theoryStatus', ...
    ['Elliptical first-order TLmean calibration obtained by multiplying ' ...
    'the feasible general-F TDmean variance by an empirical radial coefficient.']);

qrad = d2cc(:);
qrad = qrad(isfinite(qrad) & qrad >= 0);
if isempty(qrad)
    ellipticalTL.reason = 'No finite complete-case radial distances are available.';
else
    meanQ = mean(qrad);
    if ~isfinite(meanQ) || meanQ <= 0
        ellipticalTL.reason = 'The mean complete-case radial distance is not positive.';
    else
        kappaF = mean(qrad./(qrad+eps0))/meanQ;
        ellipticalTL.radialMeanQ = meanQ;
        ellipticalTL.kappa = kappaF;
        ellipticalTL.lambda = kappaF;
        ellipticalTL.coefficient = kappaF;
        if ~isfinite(kappaF) || kappaF <= 0
            ellipticalTL.reason = 'The empirical elliptical radial coefficient is not finite and positive.';
        elseif ~generalF.available
            ellipticalTL.reason = generalF.reason;
        else
            ellipticalTL.available = true;
            ellipticalTL.reason = '';
            ellipticalTL.sigma2 = kappaF^2*generalF.sigma2;
            ellipticalTL.se = sqrt(ellipticalTL.sigma2/n);
            ellipticalTL.z = Tobs(2)/ellipticalTL.se;
            ellipticalTL.pvalue = erfc(abs(ellipticalTL.z)/sqrt(2));
        end
    end
end

% ---------------------------------------------------------------------
% Gaussian-information specialization retained as diagnostic.
% ---------------------------------------------------------------------
gaussian = struct;
gaussian.available = false;
gaussian.reason = '';
gaussian.VA = NaN;
gaussian.IA = JA;
gaussian.aSigma = aSigma;
gaussian.kappa = NaN;
gaussian.TDmean = struct('available',false,'calibration','gaussianInfo', ...
    'coefficient',1,'sigma2',NaN,'se',NaN,'z',NaN,'pvalue',NaN);
gaussian.TLmean = struct('available',false,'calibration','gaussianInfo', ...
    'kappa',NaN,'lambda',NaN,'coefficient',NaN,'sigma2',NaN,'se',NaN, ...
    'z',NaN,'pvalue',NaN);
gaussian.theoryStatus = ['Gaussian-information specialization; retained as ' ...
    'a diagnostic benchmark and not generally valid under non-Gaussian F.'];

if ~isfinite(rcondJA) || rcondJA <= 1e-12
    gaussian.reason = ['Observed-data Gaussian covariance information is ' ...
        'numerically singular.'];
else
    VA = aSigma'*(JA\aSigma);
    sigmaD2G = 2*p/qhat - VA;
    tolVarG = 1e-10*max(1,2*p/qhat);

    if sigmaD2G < 0 && sigmaD2G >= -tolVarG
        sigmaD2G = 0;
    end

    gaussian.VA = VA;

    if ~isfinite(sigmaD2G) || sigmaD2G < 0
        gaussian.reason = ['Gaussian-information variance is not finite and ' ...
            'nonnegative.'];
    elseif sigmaD2G <= tolVarG
        gaussian.reason = 'Gaussian-information variance is numerically degenerate.';
    else
        % Stabilized log-ratio coefficient under Q~chi2_p.
        kappa = integral(@(q) local_kappa_integrand(q,p,eps0),0,Inf, ...
            'RelTol',1e-10,'AbsTol',1e-12)/p;

        gaussian.available = true;
        gaussian.reason = '';
        gaussian.kappa = kappa;

        gaussian.TDmean.available = true;
        gaussian.TDmean.sigma2 = sigmaD2G;
        gaussian.TDmean.se = sqrt(sigmaD2G/n);
        gaussian.TDmean.z = Tobs(4)/gaussian.TDmean.se;
        gaussian.TDmean.pvalue = erfc(abs(gaussian.TDmean.z)/sqrt(2));

        gaussian.TLmean.available = true;
        gaussian.TLmean.kappa = kappa;
        gaussian.TLmean.lambda = kappa;
        gaussian.TLmean.coefficient = kappa;
        gaussian.TLmean.sigma2 = kappa^2*sigmaD2G;
        gaussian.TLmean.se = sqrt(gaussian.TLmean.sigma2/n);
        gaussian.TLmean.z = Tobs(2)/gaussian.TLmean.se;
        gaussian.TLmean.pvalue = erfc(abs(gaussian.TLmean.z)/sqrt(2));
    end
end

% Common diagnostics, irrespective of which alpha=0 branch is primary.
asympt.qhat = qhat;
asympt.aSigma = aSigma;
asympt.JA = JA;
asympt.rcondJA = rcondJA;
asympt.patternsMissing = patt;
asympt.patternCounts = patternCounts;
asympt.patternProportions = patternProportions;
asympt.generalF = generalF;
asympt.OmegaDelta = generalF.OmegaDelta;
asympt.discrepancyIFMeanNorm = generalF.deltaMeanNorm;
asympt.omnibusVarianceIdentityResidual = generalF.varianceIdentityResidual;
asympt.ellipticalTL = ellipticalTL;
asympt.gaussian = gaussian;
asympt.VA = gaussian.VA; % backward-compatible Gaussian diagnostic field

if preferGeneralF
    % Ordinary untrimmed MDP-U: general-F is primary for TDmean.
    asympt.mode = 'feasible general-F sandwich';
    asympt.primaryCalibration = 'generalF+ellipticalRadial';
    asympt.radialModel = 'general-F with elliptical TLmean radial scaling';
    asympt.theoryStatus = ['TDmean uses the feasible general-F sandwich under ' ...
        'MCAR and finite fourth moments. Under ellipticity, TLmean uses the ' ...
        'same estimator-discrepancy variance multiplied by an empirical ' ...
        'complete-case radial coefficient. Gaussian-information results are ' ...
        'reported separately as a benchmark.'];
    asympt.aggregation = 'ordinary';
    asympt.filterCutoff = NaN;
    asympt.kc = 1;
    asympt.gammaFilter = 1;

    if generalF.available
        asympt.available = true;
        asympt.reason = '';
        asympt.degenerate = false;
        asympt.baseSigmaD2 = generalF.sigma2;
        asympt.TDmean = struct('available',true,'calibration','generalF', ...
            'coefficient',1,'sigma2',generalF.sigma2,'se',generalF.se, ...
            'z',generalF.z,'pvalue',generalF.pvalue);
    else
        asympt.available = false;
        asympt.reason = generalF.reason;
        asympt.degenerate = generalF.degenerate;
        asympt.baseSigmaD2 = NaN;
        asympt.TDmean = struct('available',false,'calibration','generalF', ...
            'coefficient',1,'sigma2',NaN,'se',NaN,'z',NaN,'pvalue',NaN);
    end

    % Under ellipticity, TLmean is a scalar multiple of TDmean to first
    % order. Use the empirical radial coefficient with the general-F
    % estimator-discrepancy variance; keep the Gaussian calculation separate.
    asympt.TLmean = ellipticalTL;
    asympt.baseKappa = ellipticalTL.kappa;
    asympt.lambda = ellipticalTL.kappa;
else
    % A general-F fixed-cutoff theorem is not currently available. Preserve
    % the previous Gaussian benchmark as the primary analytical reference;
    % the caller applies the fixed-cutoff Tallis coefficients afterwards.
    asympt.mode = 'Gaussian information benchmark';
    asympt.primaryCalibration = 'gaussianInfo';
    asympt.radialModel = 'Gaussian';
    asympt.theoryStatus = ['For alpha=0 with fixed aggregation the analytical ' ...
        'benchmark remains Gaussian; the untrimmed general-F TDmean result is ' ...
        'returned in asympt.generalF but is not applied to the fixed cutoff.'];

    if gaussian.available
        asympt.available = true;
        asympt.reason = '';
        asympt.degenerate = false;
        asympt.TDmean = gaussian.TDmean;
        asympt.TLmean = gaussian.TLmean;
        asympt.baseSigmaD2 = gaussian.TDmean.sigma2;
        asympt.baseKappa = gaussian.kappa;
    else
        asympt.available = false;
        asympt.reason = gaussian.reason;
        asympt.degenerate = false;
    end
end
end

% -------------------------------------------------------------------------
function Calibration = local_calibration_labels(alpha,aggregation, ...
    consistencyfactor,coupledtrim)
%local_calibration_labels Labels used in the human-readable results table.

Calibration = {'not available';'not available';'not available';'not available'};

if alpha == 0
    if strcmp(aggregation,'ordinary')
        Calibration{2} = 'elliptical radial sandwich';
        Calibration{4} = 'general-F sandwich';
    elseif strcmp(aggregation,'fixed')
        Calibration{2} = 'Gaussian fixed-cutoff benchmark';
        Calibration{4} = 'Gaussian fixed-cutoff benchmark';
    end
else
    if strcmp(consistencyfactor,'adaptive')
        if strcmp(aggregation,'inlier')
            Calibration{2} = 'adaptive generic retained-moment sandwich';
            Calibration{4} = 'adaptive generic retained-moment sandwich';
        elseif strcmp(aggregation,'fixed')
            Calibration{2} = 'adaptive elliptical fixed-cutoff';
            Calibration{4} = 'adaptive elliptical fixed-cutoff';
        else
            Calibration{2} = 'adaptive elliptical sandwich';
            Calibration{4} = 'adaptive elliptical sandwich';
        end
    elseif strcmp(consistencyfactor,'pattern')
        if coupledtrim
            Calibration{2} = ['Gaussian-pattern generic retained-moment TEM sandwich ' ...
                '(frozen eligibility)'];
            Calibration{4} = ['Gaussian-pattern generic retained-moment TEM sandwich ' ...
                '(frozen eligibility)'];
        elseif strcmp(aggregation,'inlier')
            Calibration{2} = 'Gaussian-pattern generic retained-moment sandwich';
            Calibration{4} = 'Gaussian-pattern generic retained-moment sandwich';
        elseif strcmp(aggregation,'fixed')
            Calibration{2} = 'Gaussian fixed-cutoff sandwich';
            Calibration{4} = 'Gaussian fixed-cutoff sandwich';
        else
            Calibration{2} = 'Gaussian pattern sandwich';
            Calibration{4} = 'Gaussian pattern sandwich';
        end
    end
end
end

% -------------------------------------------------------------------------
function y = local_kappa_integrand(q,p,eps0)
%local_kappa_integrand Integrand for E{Q/(Q+eps0)}, Q~chi2_p.

ratio = q./(q+eps0);
ratio(q == 0) = 0;
ratio(isinf(q)) = 1;
y = ratio.*chi2pdf(q,p);
y(~isfinite(y)) = 0;
end

% -------------------------------------------------------------------------
function asymptStar = local_bootstrap_omnibus_asymptotic( ...
    Ystar,maskMiss,completeIdx,outFitStar,alpha,method,robustClass, ...
    robustEff,robustBonflev,Tstar,eps0,outRobStar,coupledtrim, ...
    forcedZeroStar,consistencyfactor,muCCStar,SigCCStar,adaptivepool, ...
    adaptiveminref,d2ccStar,d2allccStar,muHatStar,SigHatStar)
%local_bootstrap_omnibus_asymptotic Replicate-specific discrepancy sandwich.
%
% This helper reconstructs the same end-to-end scatter-discrepancy covariance
% estimator used for the observed omnibus statistic, but evaluated entirely in
% one bootstrap replicate. Final scalar aggregation does not enter because the
% omnibus statistic acts on the fitted scatter-estimator discrepancy itself.

n = size(Ystar,1);
nCompleteStar = sum(completeIdx);

if coupledtrim
    % Frozen-eligibility working sandwich, matching the observed coupled branch.
    asymptStar = local_tem_asymptotic(Ystar,maskMiss,completeIdx,outFitStar, ...
        alpha,method,robustClass,robustEff,robustBonflev,Tstar,eps0, ...
        outRobStar,~forcedZeroStar);
elseif alpha == 0
    % The omnibus discrepancy covariance is the general-F covariance even when
    % the final scalar aggregation is not ordinary.
    TDordinaryStar = mean(d2allccStar-d2ccStar);
    asymptStar = local_classical_asymptotic(Ystar,maskMiss,completeIdx, ...
        muHatStar,SigHatStar,nCompleteStar,Tstar,TDordinaryStar,eps0, ...
        d2ccStar,true);
elseif strcmp(consistencyfactor,'pattern')
    asymptStar = local_tem_asymptotic(Ystar,maskMiss,completeIdx,outFitStar, ...
        alpha,method,robustClass,robustEff,robustBonflev,Tstar,eps0, ...
        outRobStar,[]);
elseif strcmp(consistencyfactor,'adaptive')
    asymptStar = local_tem_adaptive_asymptotic(Ystar,maskMiss,completeIdx, ...
        outFitStar,alpha,method,robustClass,robustEff,robustBonflev, ...
        Tstar,eps0,muCCStar,SigCCStar,adaptivepool,adaptiveminref,outRobStar);
else
    asymptStar = local_unavailable_asymptotic(nCompleteStar/n);
    asymptStar.reason = ['Bootstrap omnibus covariance is unavailable for ' ...
        'the selected TEM consistency-factor option.'];
end
end

% -------------------------------------------------------------------------
function omnibus = local_omnibus_test(asympt,SigA,SigC,n,alpha, ...
    consistencyfactor,coupledtrim,rankExponent,nEff,nEffSource)
%local_omnibus_test Rank-selected Hausman scatter-discrepancy diagnostic.
%
% Let delta=vech(SigA-SigC) and let OmegaDelta estimate the asymptotic
% covariance of sqrt(n)*delta. If lambda_1 >= ... >= lambda_s >= 0 are the
% eigenvalues of OmegaDelta, define
%
%   eta_n = nEff^(-rankExponent),
%   tau_n = eta_n*lambda_1,
%
% and retain eigenvalue j when lambda_j/lambda_1 > eta_n. The statistic is
%
%   H = n*delta'*OmegaDelta_eta^+*delta,
%
% where OmegaDelta_eta^+ inverts only the statistically retained spectral
% subspace. Its analytical null reference is chi-square with degrees of
% freedom equal to the selected rank. A machine-precision tolerance is used
% only to diagnose numerical non-positive-semidefiniteness; it does not
% determine the statistical rank.

p = size(SigA,1);
s = p*(p+1)/2;
omnibus = struct('available',false,'reason','', ...
    'stat',NaN,'df',NaN,'pvalue',NaN,'rank',0,'dimension',s, ...
    'delta',NaN(s,1),'OmegaDelta',NaN(s,s),'OmegaDeltaPlus',NaN(s,s), ...
    'eigenvalues',NaN(s,1),'relativeEigenvalues',NaN(s,1), ...
    'largestEigenvalue',NaN,'numericalTolerance',NaN, ...
    'rankTolerance',NaN,'rankThreshold',NaN,'relativeRankThreshold',NaN, ...
    'rankExponent',rankExponent,'effectiveSampleSize',nEff, ...
    'effectiveSampleSizeSource',nEffSource,'rhoLastKept',NaN, ...
    'rhoFirstDiscarded',NaN,'spectralGapRatio',NaN, ...
    'lowerThresholdMargin',NaN,'upperThresholdMargin',NaN, ...
    'rankSensitivity',table(), ...
    'nullSpaceResidual',NaN,'calibration','not available', ...
    'theoryStatus',['Rank-selected Hausman-type quadratic estimator-discrepancy ' ...
    'test for the complete-case versus all-available scatter estimators.']);

if ~isequal(size(SigA),size(SigC)) || size(SigA,2)~=p || ...
        any(~isfinite(SigA(:))) || any(~isfinite(SigC(:)))
    omnibus.reason = 'The fitted scatter matrices are invalid.';
    return
end

SA = (SigA+SigA')/2;
SC = (SigC+SigC')/2;
deltaHat = local_vech(SA-SC);
omnibus.delta = deltaHat;

if ~isstruct(asympt) || ~isfield(asympt,'OmegaDelta') || ...
        ~isnumeric(asympt.OmegaDelta) || ...
        ~isequal(size(asympt.OmegaDelta),[s s]) || ...
        any(~isfinite(asympt.OmegaDelta(:)))
    omnibus.reason = ['The full scatter-discrepancy influence covariance is ' ...
        'not available for this analytical configuration.'];
    return
end

Omega = (asympt.OmegaDelta+asympt.OmegaDelta')/2;
omnibus.OmegaDelta = Omega;

% Eigenvalues and a purely numerical tolerance. The numerical tolerance is
% deliberately separated from the statistical spectral threshold below.
[V,Dmat] = eig(Omega);
D = real(diag(Dmat));
V = real(V);
maxAbsEig = max(abs(D));
if isempty(maxAbsEig) || ~isfinite(maxAbsEig) || maxAbsEig <= 0
    omnibus.reason = 'The estimated discrepancy covariance is numerically zero.';
    return
end
numericalTol = max(s,1)*eps(maxAbsEig);
omnibus.numericalTolerance = numericalTol;
% Backward-compatible field name: it now stores only the numerical tolerance,
% not the statistical rank threshold.
omnibus.rankTolerance = numericalTol;

% Tiny negative eigenvalues can arise from roundoff. A materially negative
% eigenvalue signals an invalid estimated sandwich matrix.
if any(D < -100*numericalTol)
    omnibus.eigenvalues = sort(D,'descend');
    omnibus.reason = ['The estimated discrepancy covariance is not positive ' ...
        'semidefinite within numerical tolerance.'];
    return
end
D(D < 0) = 0;

% Sort the eigensystem in descending order so that boundary ratios and the
% sensitivity table have a deterministic interpretation.
[D,ord] = sort(D,'descend');
V = V(:,ord);
lambda1 = D(1);
if ~isfinite(lambda1) || lambda1 <= numericalTol
    omnibus.eigenvalues = D;
    omnibus.reason = 'The estimated discrepancy covariance has numerical rank zero.';
    return
end

rhoEig = D/lambda1;
eta = nEff^(-rankExponent);
rankThreshold = eta*lambda1;
keep = rhoEig > eta;
r = sum(keep);

omnibus.eigenvalues = D;
omnibus.relativeEigenvalues = rhoEig;
omnibus.largestEigenvalue = lambda1;
omnibus.relativeRankThreshold = eta;
omnibus.rankThreshold = rankThreshold;
omnibus.rank = r;
omnibus.df = r;

if r >= 1
    omnibus.rhoLastKept = rhoEig(r);
    omnibus.lowerThresholdMargin = omnibus.rhoLastKept/eta;
end
if r < s
    omnibus.rhoFirstDiscarded = rhoEig(r+1);
    if omnibus.rhoFirstDiscarded > 0
        omnibus.upperThresholdMargin = eta/omnibus.rhoFirstDiscarded;
        if r >= 1
            omnibus.spectralGapRatio = ...
                omnibus.rhoLastKept/omnibus.rhoFirstDiscarded;
        end
    elseif r >= 1 && omnibus.rhoFirstDiscarded == 0
        omnibus.upperThresholdMargin = Inf;
        omnibus.spectralGapRatio = Inf;
    end
end

% Multiplicative rank-separation diagnostics. Values comfortably above one
% indicate that the last retained eigenvalue lies above the threshold and the
% first discarded eigenvalue lies below it. spectralGapRatio compares the two
% boundary eigenvalue ratios directly. These are descriptive finite-sample
% diagnostics; they do not alter the rank-selection rule.

% Finite-sample sensitivity requested in the paper: eta_n=nEff^(-a) for
% a in {1/4,1/3,2/5}. Larger a gives a smaller threshold and hence a weakly
% larger selected rank.
expSens = [1/4;1/3;2/5];
etaSens = nEff.^(-expSens);
thrSens = etaSens*lambda1;
rankSens = zeros(numel(expSens),1);
for jj = 1:numel(expSens)
    rankSens(jj) = sum(rhoEig > etaSens(jj));
end
omnibus.rankSensitivity = table(expSens,etaSens,thrSens,rankSens, ...
    'VariableNames',{'Exponent','Eta','Threshold','Rank'});

if r == 0
    omnibus.reason = ['The statistical spectral threshold retains no ' ...
        'scatter-discrepancy direction.'];
    return
end

Vr = V(:,keep);
dr = D(keep);
OmegaPlus = Vr*diag(1./dr)*Vr';
OmegaPlus = (OmegaPlus+OmegaPlus')/2;
omnibus.OmegaDeltaPlus = OmegaPlus;

H = n*(deltaHat'*OmegaPlus*deltaHat);
if H < 0 && H >= -100*eps(max(1,abs(H)))
    H = 0;
end
if ~isfinite(H) || H < 0
    omnibus.reason = 'The omnibus quadratic statistic is not finite and nonnegative.';
    return
end

% Component of the observed discrepancy outside the statistically selected
% covariance range. It is not included in the thresholded quadratic form.
projDelta = Vr*(Vr'*deltaHat);
omnibus.nullSpaceResidual = norm(deltaHat-projDelta)/max(norm(deltaHat),eps);

omnibus.stat = H;
omnibus.pvalue = gammainc(H/2,r/2,'upper');
omnibus.available = isfinite(omnibus.pvalue);
if ~omnibus.available
    omnibus.reason = 'Unable to evaluate the chi-square upper-tail probability.';
    return
end

if alpha == 0
    omnibus.calibration = 'rank-selected general-F scatter-discrepancy sandwich';
elseif strcmp(consistencyfactor,'adaptive')
    omnibus.calibration = 'rank-selected adaptive elliptical scatter-discrepancy sandwich';
elseif strcmp(consistencyfactor,'pattern') && coupledtrim
    omnibus.calibration = ['rank-selected Gaussian pattern scatter-discrepancy ' ...
        'sandwich (frozen eligibility)'];
elseif strcmp(consistencyfactor,'pattern')
    omnibus.calibration = 'rank-selected Gaussian pattern scatter-discrepancy sandwich';
else
    omnibus.calibration = 'rank-selected scatter-discrepancy sandwich';
end
omnibus.reason = '';
end

% -------------------------------------------------------------------------
function asympt = local_unavailable_asymptotic(qhat)
%local_unavailable_asymptotic Initialize analytical benchmark output.

asympt = struct;
asympt.available = false;
asympt.reason = 'Analytical benchmark is not available for this configuration.';
asympt.jointReason = '';
asympt.qhat = qhat;
asympt.VA = NaN;
asympt.degenerate = false;
asympt.aggregation = '';
asympt.filterCutoff = NaN;
asympt.kc = NaN;
asympt.lambda = NaN;
asympt.gammaFilter = NaN;
asympt.baseSigmaD2 = NaN;
asympt.baseKappa = NaN;
asympt.OmegaDelta = NaN;
asympt.OmegaDeltaMu = NaN;
asympt.OmegaDeltaMuSigma = NaN;
asympt.OmegaJoint = NaN;
asympt.discrepancyIFMeanNorm = NaN;
asympt.locationDiscrepancyIFMeanNorm = NaN;
asympt.omnibusVarianceIdentityResidual = NaN;
[asympt.inlierScaling,~] = local_stable_inlier_scaling();
asympt.TDmean = struct('coefficient',NaN,'sigma2',NaN,'se',NaN, ...
    'z',NaN,'pvalue',NaN);
asympt.TLmean = struct('kappa',NaN,'lambda',NaN,'coefficient',NaN, ...
    'sigma2',NaN,'se',NaN,'z',NaN,'pvalue',NaN);
end

% -------------------------------------------------------------------------
function D = local_duplication_matrix(p)
%local_duplication_matrix Duplication matrix for column-wise vech.
%
% D satisfies vec(H)=D*vech(H) for every symmetric p-by-p matrix H, where
% vech stacks the lower triangular part column by column.

s = p*(p+1)/2;
D = zeros(p*p,s);
k = 0;
for j = 1:p
    for i = j:p
        k = k + 1;
        E = zeros(p,p);
        E(i,j) = 1;
        if i ~= j
            E(j,i) = 1;
        end
        D(:,k) = E(:);
    end
end
end

% -------------------------------------------------------------------------
function U = local_random_spherical_directions(n,p)
%local_random_spherical_directions iid uniform directions on S^(p-1).
%
% If Z~N_p(0,I), then Z/||Z|| is uniform on the unit sphere and is
% independent of ||Z||. The Gaussian draw is therefore used only to obtain
% rotationally invariant directions; it does not impose a Gaussian radial
% distribution on the empirical-radial bootstrap sample.

U = randn(n,p);
normU = sqrt(sum(U.^2,2));
bad = ~isfinite(normU) | normU <= 0;
while any(bad)
    U(bad,:) = randn(sum(bad),p);
    normU(bad) = sqrt(sum(U(bad,:).^2,2));
    bad = ~isfinite(normU) | normU <= 0;
end
U = U./normU;
end

% -------------------------------------------------------------------------
function Sspd = local_make_spd(S)
% Make covariance matrix symmetric positive definite if needed.

S = (S + S')/2;
[~,flag] = chol(S);

if flag == 0
    Sspd = S;
    return
end

lam = 1e-8 * trace(S) / size(S,1);
if lam <= 0 || ~isfinite(lam)
    lam = 1e-8;
end

I = eye(size(S));
for k = 1:8
    Stmp = S + lam*I;
    [~,flag] = chol(Stmp);
    if flag == 0
        Sspd = Stmp;
        return
    end
    lam = 10*lam;
end

error('FSDA:mdMDPtest:NonSPD', ...
    'Unable to regularize covariance matrix to positive definiteness.');
end


% -------------------------------------------------------------------------
function a = local_tem_inv_adjust(c,pg,p,n,method)
%local_tem_inv_adjust Inverse parameter-free adjusted-distance map.
%
% The analytical alpha>0 sandwich is built on the raw p_g-dimensional
% squared Mahalanobis radius Q_g. TEM trims on an adjusted full-dimensional
% distance m_g(Q_g), so this helper returns a_g(c)=m_g^{-1}(c).
%
% For the supported maps the transformation depends on the missingness
% pattern only through p_g (and, for betaMap, n), not through Sigma. This is
% why the existing scatter Jacobian remains valid after replacing the pri
% cutoff by the appropriate a_g(c). detMap is intentionally excluded because
% its transformation depends explicitly on the full and observed-block
% determinants of Sigma.

method = string(method);
switch method
    case "pri"
        a = c-(p-pg);

    case "expScale"
        a = c*pg/p;

    case "zMap"
        a = pg+sqrt(pg/p)*(c-p);

    case "chiMap"
        u = chi2cdf(c,p);
        u = min(max(u,eps),1-eps);
        a = chi2inv(u,pg);

    case "betaMap"
        cn = (n-1)^2/n;
        if n <= p+1 || n <= pg+1
            a = c;
            return
        end
        u = min(max(c/cn,0),1-eps);
        al = betacdf(u,p/2,(n-p-1)/2);
        al = min(max(al,eps),1-eps);
        a = cn*betainv(al,pg/2,(n-pg-1)/2);

    otherwise
        a = NaN;
end
end

% -------------------------------------------------------------------------
function d2adj = local_adjusted_distances(Y,mu,Sigma,method)
%local_adjusted_distances Final adjusted distances for analytical fallbacks.
%
% This helper mirrors the parameter-free distance calculations in mdTEM.
% detMap is deliberately not handled because the analytical branches return
% before reaching this fallback.

p = size(Y,2);
[d2part,pobs] = mdPartialMD(Y,mu,Sigma);
d2adj = mdPartialMD2full(d2part,p,pobs,'method',method);
d2adj = d2adj(:);
end

% -------------------------------------------------------------------------
function out = local_coupled_tem(Y,forcedZero,alpha,method,tol)
%local_coupled_tem Pattern-corrected TEM with fixed forced-zero rows.
%
% This local sensitivity implementation mirrors the pattern-corrected mdTEM
% concentration step. At each iteration it first constructs TEM's ordinary
% indicator wT from the h=floor(n*(1-alpha)) smallest adjusted distances and
% then applies the fixed eligibility constraint
%
%   w = wT .* (~forcedZero).
%
% Consequently a complete case rejected by the robust complete-case fit can
% never contribute to the TEM E/M step. The internal threshold is still the
% h-th ordinary TEM adjusted-distance threshold. Because the extra exclusion
% is not generated solely by that radial threshold, the usual pattern-wise
% Tallis correction is therefore a working correction in this sensitivity
% mode. mdMDPtest reports a frozen-eligibility analytical approximation in
% which this externally constructed eligibility mask is treated as fixed.

[n,p] = size(Y);
forcedZero = logical(forcedZero(:));
if numel(forcedZero) ~= n
    error('FSDA:mdMDPtest:WrongCoupledMask', ...
        'The coupled TEM mask must have one entry for every row of Y.');
end

maxiter = 100;
method = string(method);

% Same initialization used by mdTEM.
mus = mean(Y,1,"omitmissing")';
X0 = Y;
for j = 1:p
    miss = isnan(X0(:,j));
    X0(miss,j) = mus(j);
end
sigs = cov(X0,1);

keep_count = max(0,floor(n*(1-alpha)));
nanY = isnan(Y);
w = zeros(n,1);
wT = zeros(n,1);
kinfo = [];
kfactor = 1;
dif = Inf;
iter = 0;

while dif > tol && iter < maxiter
    iter = iter + 1;
    mus_old = mus;
    sigs_old = sigs;

    if method == "impMD"
        Yimp = mdImputeCondMean(Y,mus,sigs);
        d2_adj = mahalFS(Yimp,mus',sigs);
        poss = sum(~nanY,2);
    elseif method == "detMap"
        [d2,poss] = mdPartialMD(Y,mus,sigs);
        d2_adj = mdPartialMD2full(d2,p,poss,'method',method, ...
            'Y',Y,'Sigma',sigs);
    else
        [d2,poss] = mdPartialMD(Y,mus,sigs);
        d2_adj = mdPartialMD2full(d2,p,poss,'method',method);
    end

    nanMask = isnan(d2_adj);
    [~,idxSorted] = sort(d2_adj,'ascend','MissingPlacement','last');
    keepIdx = idxSorted(1:min(keep_count,sum(~nanMask)));
    if isempty(keepIdx)
        error('FSDA:mdMDPtest:NoTEMRows', ...
            'No finite row is available for the coupled TEM concentration step.');
    end

    % Ordinary TEM concentration indicator and coupled product weight.
    wT = zeros(n,1);
    wT(keepIdx) = 1;
    w = wT;
    w(forcedZero) = 0;
    if ~any(w)
        error('FSDA:mdMDPtest:NoTEMRows', ...
            'The coupled TEM constraint removes all retained rows.');
    end

    % Keep the ordinary TEM threshold. The additional fixed exclusion is
    % represented only through w, exactly as in the product definition.
    cthr = max(d2_adj(keepIdx));

    kinfo = local_coupled_pattern_factors(nanY,w,poss,cthr,p,n,sigs,method);
    [mus,sigs,kfactor] = local_coupled_corrected_step(Y,nanY,w,mus, ...
        sigs,kinfo);

    muDiff = max(abs(mus(:)-mus_old(:)));
    sigmaDiff = max(abs(sigs(:)-sigs_old(:)));
    dif = max(muDiff,sigmaDiff);
end

out = struct;
out.loc = mus;
out.cov = sigs;
out.iter = iter;
out.weights = w;
out.internalWeights = wT;
out.forcedZero = forcedZero;
out.kfactor = kfactor;
out.kinfo = kinfo;
out.cthr = cthr;
out.adjustedD2 = d2_adj;
end

% -------------------------------------------------------------------------
function kinfo = local_coupled_pattern_factors(nanY,w,poss,cthr,p,n,Sigma,method)
%local_coupled_pattern_factors Pattern-wise Tallis factors for coupled TEM.

keep = w > 0;
patt = unique(nanY(keep,:),'rows');
G = size(patt,1);
pobs = zeros(G,1);
nkept = zeros(G,1);
athr = zeros(G,1);
gammag = zeros(G,1);
kg = nan(G,1);

for g = 1:G
    pg = patt(g,:);
    rows = keep & all(nanY == pg,2);
    obs = ~pg;
    pobs(g) = sum(obs);
    nkept(g) = sum(rows);

    if pobs(g) == 0
        athr(g) = 0;
        gammag(g) = 0;
        kg(g) = NaN;
        continue
    end

    athr(g) = local_coupled_inv_adjust(cthr,pobs(g),p,n,Sigma,obs,method);
    if athr(g) <= 0 || ~isfinite(athr(g))
        gammag(g) = 0;
        kg(g) = NaN;
    else
        gammag(g) = chi2cdf(athr(g),pobs(g));
        if gammag(g) <= 0
            kg(g) = NaN;
        else
            kg(g) = chi2cdf(athr(g),pobs(g)+2)/gammag(g);
        end
    end
end

kinfo = table(pobs,nkept,athr,gammag,kg);
kbar = local_coupled_weighted_factor(kinfo);
bad = ~isfinite(kinfo.kg) | kinfo.kg <= 0 | kinfo.kg > 1 | ...
    kinfo.nkept < max(2,kinfo.pobs);
kinfo.kg(bad) = kbar;
end

% -------------------------------------------------------------------------
function kbar = local_coupled_weighted_factor(kinfo)
%local_coupled_weighted_factor Information-weighted pattern factor.

ok = isfinite(kinfo.kg) & kinfo.kg > 0 & kinfo.kg <= 1 & kinfo.nkept > 0;
if ~any(ok)
    kbar = 1;
    return
end
wgt = kinfo.nkept(ok).*kinfo.pobs(ok);
kbar = sum(wgt.*kinfo.kg(ok))/sum(wgt);
if ~isfinite(kbar) || kbar <= 0
    kbar = 1;
end
end

% -------------------------------------------------------------------------
function a = local_coupled_inv_adjust(c,pg,p,n,Sigma,obs,method)
%local_coupled_inv_adjust Inverse adjusted-distance map for one pattern.

method = string(method);
switch method
    case "pri"
        a = c-(p-pg);
    case "expScale"
        a = c*pg/p;
    case "zMap"
        a = pg+sqrt(pg/p)*(c-p);
    case "detMap"
        gfull = exp(local_coupled_logdet_spd(Sigma)/p);
        gobs = exp(local_coupled_logdet_spd(Sigma(obs,obs))/pg);
        a = c*(pg/p)*(gobs/gfull);
    case "chiMap"
        u = chi2cdf(c,p);
        u = min(max(u,eps),1-eps);
        a = chi2inv(u,pg);
    case "betaMap"
        cn = (n-1)^2/n;
        if n <= p+1 || n <= pg+1
            a = c;
            return
        end
        u = min(max(c/cn,0),1-eps);
        al = betacdf(u,p/2,(n-p-1)/2);
        al = min(max(al,eps),1-eps);
        a = cn*betainv(al,pg/2,(n-pg-1)/2);
    case "impMD"
        a = c;
    otherwise
        a = c;
end
end

% -------------------------------------------------------------------------
function [musNew,sigsNew,kbar] = local_coupled_corrected_step(Y,nanY,w, ...
    mus,sigs,kinfo)
%local_coupled_corrected_step Pattern-wise corrected E/M update.

p = size(Y,2);
mus = mus(:);
keep = w > 0;
S = zeros(p,p);
m1 = zeros(p,1);
h = 0;
patt = unique(nanY(keep,:),'rows');

for g = 1:size(patt,1)
    pg = patt(g,:);
    rows = keep & all(nanY == pg,2);
    hg = sum(rows);
    if hg == 0
        continue
    end

    o = find(~pg);
    m = find(pg);
    if isempty(o)
        continue
    end

    idx = find(kinfo.pobs == numel(o) & kinfo.nkept == hg,1);
    if isempty(idx)
        kg = 1;
    else
        kg = kinfo.kg(idx);
    end
    if ~isfinite(kg) || kg <= 0
        kg = 1;
    end

    Z = Y(rows,o)-mus(o)';
    Szz = (Z'*Z)/kg;
    sZ = sum(Z,1)';

    S(o,o) = S(o,o)+Szz;
    m1(o) = m1(o)+sZ;

    if ~isempty(m)
        Soo = sigs(o,o);
        Ag = sigs(m,o)/Soo;
        Cg = sigs(m,m)-Ag*sigs(o,m);
        Cg = (Cg+Cg')/2;

        SzzAg = Szz*Ag';
        S(o,m) = S(o,m)+SzzAg;
        S(m,o) = S(m,o)+SzzAg';
        S(m,m) = S(m,m)+Ag*SzzAg+hg*Cg;
        m1(m) = m1(m)+Ag*sZ;
    end
    h = h+hg;
end

if h == 0
    musNew = mus;
    sigsNew = sigs;
    kbar = 1;
    return
end

m1 = m1/h;
S = S/h;
musNew = mus+m1;
sigsNew = S-(m1*m1');
sigsNew = (sigsNew+sigsNew')/2;
kbar = local_coupled_weighted_factor(kinfo);
end

% -------------------------------------------------------------------------
function val = local_coupled_logdet_spd(S)
%local_coupled_logdet_spd Stable log determinant for a positive scatter matrix.

S = (S+S')/2;
[R,flag] = chol(S);
if flag == 0
    val = 2*sum(log(diag(R)));
else
    val = log(max(det(S),realmin));
end
end

% -------------------------------------------------------------------------
function so = mdMDPsecondOrder(p,alpha)
%mdMDPsecondOrder Complete-Gaussian adaptive-TEM second-order benchmark.
%
%  so = mdMDPsecondOrder(p,alpha) returns the final O(n^{-1}) centering
%  constants for the complete-data Gaussian benchmark with the adaptive
%  radial TEM correction.  The decomposition is
%
%        cTotal = cLin + cQuad + cMu,
%
%  where cMu is the location contribution, cQuad is the inverse-scatter
%  curvature contribution and cLin is the moving-boundary/adaptive-feedback
%  contribution.
%
%  IMPORTANT.  This local function is a theoretical benchmark, not a
%  general second-order correction for data with missing values.  In
%  particular, the present closed form for cLin is derived for the complete
%  Gaussian adaptive benchmark.  mdMDPtest therefore reports these constants
%  in out.asympt.secondOrderBenchmark but does not use them to modify p-values.
%
%  Input arguments:
%
%    p     : Dimension. Positive integer >= 2.
%    alpha : Trimming fraction in [0,0.5].
%
%  Output fields:
%
%    p, alpha, gamma : Dimension, trimming fraction and retained fraction.
%    a               : chi-square_p gamma quantile.
%    f               : chi-square_p density at a.
%    beta            : F_{p+2}(a).
%    js              : F_{p+4}(a), trace-free/shape Jacobian coefficient.
%    jt              : Scalar trace Jacobian coefficient,
%                       ((p+2)*gamma*js-p*beta^2)/(2*gamma).
%    kappa           : Population Gaussian Tallis factor beta/gamma.
%    lambda          : Adaptive scalar feedback coefficient, 1-jt/beta.
%    feedback1       : 1/(1-lambda).
%    feedback2       : 1/(1-lambda)^2.
%    cMu             : Location contribution.
%    cQuad           : Inverse-scatter curvature contribution.
%    cLin            : Moving-boundary/adaptive-feedback contribution.
%    cTotal          : Total complete-Gaussian adaptive centering constant.
%    cLinTerms       : 1 x 3 vector containing the three terms of cLin.
%    scope           : Text identifying the theoretical scope.
%    usedForCalibration : false.  The constants are diagnostic only.
%
%  The chi-square recurrences
%
%      gamma-beta = 2*a*f/p,
%      beta-js    = 2*a^2*f/[p(p+2)]
%
%  are used for numerically stable evaluation of the small differences in
%  the cLin formula.  The algebraically equivalent closed form is
%
%    cLin = p*(1/jt-1)
%         + (js-jt)/(2*jt)*cQuad
%         + (beta-jt)*((p-1)*(p+2)*jt*(beta-js)-2*beta) ...
%           /(2*jt^2*(gamma-beta)).
%
%  For p=5 and alpha=0.25 the function returns approximately
%
%      cMu    =   4.413807067927
%      cQuad  =  58.426913409290
%      cLin   =  41.759878370847
%      cTotal = 104.600598848064.
%

if ~isscalar(p) || ~isnumeric(p) || ~isfinite(p) || p < 2 || p ~= floor(p)
    error('FSDA:mdMDPtest:WrongSecondOrderP', ...
        'The second-order benchmark requires an integer p >= 2.');
end
if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || ...
        alpha < 0 || alpha > 0.5
    error('FSDA:mdMDPtest:WrongSecondOrderAlpha', ...
        'The second-order benchmark requires alpha in [0,0.5].');
end

gamma = 1-alpha;

% At alpha=0 the complete-data adaptive TEM fit coincides with the
% untrimmed reference fit.  Return the continuous no-trimming limit directly
% and avoid the removable 0/0 terms in the positive-trimming formula.
if alpha == 0
    a = Inf;
    f = 0;
    beta = 1;
    js = 1;
    jt = 1;
    kappa = 1;
    lambda = 0;
    feedback1 = 1;
    feedback2 = 1;
    cMu = 0;
    cQuad = 0;
    cLinTerms = [0 0 0];
    cLin = 0;
    cTotal = 0;
else
    a = chi2inv(gamma,p);
    f = chi2pdf(a,p);
    beta = chi2cdf(a,p+2);
    js = chi2cdf(a,p+4);

    % Scalar trace Jacobian.  This is algebraically equivalent to
    % beta+a*f*(beta/gamma-a/p), but the CDF form mirrors the final theory.
    jt = ((p+2)*gamma*js-p*beta^2)/(2*gamma);
    kappa = beta/gamma;

    if ~isfinite(jt) || jt <= 0 || ~isfinite(beta) || beta <= 0
        error('FSDA:mdMDPtest:SecondOrderUnstable', ...
            ['The complete-Gaussian adaptive second-order benchmark is ' ...
            'outside the regular stable regime (nonpositive jt or beta).']);
    end

    lambda = 1-jt/beta;
    if ~isfinite(lambda) || lambda >= 1
        error('FSDA:mdMDPtest:SecondOrderNoncontractive', ...
            ['The complete-Gaussian adaptive second-order benchmark is ' ...
            'not contractive because lambda >= 1.']);
    end

    feedback1 = 1/(1-lambda);
    feedback2 = feedback1^2;

    cMu = p*(1/beta-1);
    cQuad = (p-1)*(p+2)*(1/js-1);

    % Use the recurrence forms of gamma-beta and beta-js in the small
    % differences.  This is more stable than subtracting nearby CDF values
    % when alpha is small, while remaining exactly the same theoretically.
    gammaMinusBeta = 2*a*f/p;
    betaMinusJs = 2*a^2*f/(p*(p+2));

    if ~isfinite(gammaMinusBeta) || gammaMinusBeta <= 0
        error('FSDA:mdMDPtest:SecondOrderBoundary', ...
            'Unable to evaluate the positive-trimming boundary coefficient.');
    end

    term1 = p*(1/jt-1);
    term2 = ((js-jt)/(2*jt))*cQuad;
    term3 = ((beta-jt) * ...
        ((p-1)*(p+2)*jt*betaMinusJs-2*beta)) / ...
        (2*jt^2*gammaMinusBeta);

    cLinTerms = [term1 term2 term3];
    cLin = sum(cLinTerms);
    cTotal = cLin+cQuad+cMu;
end

so = struct;
so.p = p;
so.alpha = alpha;
so.gamma = gamma;
so.a = a;
so.f = f;
so.beta = beta;
so.js = js;
so.jt = jt;
so.kappa = kappa;
so.lambda = lambda;
so.feedback1 = feedback1;
so.feedback2 = feedback2;
so.cMu = cMu;
so.cQuad = cQuad;
so.cLin = cLin;
so.cTotal = cTotal;
so.cLinTerms = cLinTerms;
so.scope = 'complete-Gaussian adaptive-TEM second-order benchmark';
so.usedForCalibration = false;
end

%FScategory:MULT-MissingData