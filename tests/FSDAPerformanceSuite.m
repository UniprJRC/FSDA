classdef FSDAPerformanceSuite < matlab.perftest.TestCase
%FSDAPerformanceSuite Comprehensive performance tests for FSDA toolbox
%
%   Unified test suite covering all core FSDA functions across regression,
%   multivariate, and clustering domains at Small/Medium/Large data sizes.
%
%   Functions tested (13 total):
%     Regression:    FSR, FSRmdr, Sreg, MMreg, LXS (LMS), LXS (LTS)
%     Multivariate:  FSM, FSMmmd, Smult, MMmult, mcd, mahalFS
%     Clustering:    tclust, tkmeans
%
%   Run with: results = runperf('FSDAPerformanceSuite');
%
%   Authors: Anoush Najarian & Sindhuja Parimalarangan, June 2026

    properties (MethodSetupParameter)
        DataSize = struct( ...
            'Small',  struct('n', 500,  'v', 5,  'p', 5,  'k', 3), ...
            'Medium', struct('n', 2000, 'v', 10, 'p', 10, 'k', 4), ...
            'Large',  struct('n', 5000, 'v', 15, 'p', 15, 'k', 5))
    end

    properties
        Y           % Multivariate data (n x v)
        X           % Regression predictors (n x p-1)
        y           % Regression response (n x 1)
        bsb         % Initial subset for forward search
        mu          % Sample mean of Y
        Sigma       % Sample covariance of Y
    end

    methods (TestMethodSetup)
        function setupData(testCase, DataSize)
            rng(42, 'twister');
            n = DataSize.n;
            v = DataSize.v;
            p = DataSize.p;
            k = DataSize.k;

            % Multivariate data with k clusters + 10% outliers
            nPerCluster = floor(n / k);
            Ymat = zeros(n, v);
            for j = 1:k
                startIdx = (j-1)*nPerCluster + 1;
                if j == k
                    endIdx = n;
                else
                    endIdx = j*nPerCluster;
                end
                clusterSize = endIdx - startIdx + 1;
                mu_j = 5 * randn(1, v);
                Ymat(startIdx:endIdx, :) = randn(clusterSize, v) + mu_j;
            end
            nOut = round(0.1 * n);
            Ymat(1:nOut, :) = Ymat(1:nOut, :) + 15;
            testCase.Y = Ymat;

            % Regression data with 10% leverage outliers
            Xmat = randn(n, p-1);
            beta = randn(p, 1);
            yVec = [ones(n,1) Xmat] * beta + randn(n, 1);
            yVec(1:nOut) = yVec(1:nOut) + 30;
            testCase.X = Xmat;
            testCase.y = yVec;

            % Initial subset for forward search
            testCase.bsb = (1:p+1)';

            % Covariance for mahalFS
            testCase.mu = mean(Ymat);
            testCase.Sigma = cov(Ymat);
        end
    end

    methods (Test)
        %% --- Regression ---

        function testFSR(testCase, DataSize) %#ok<INUSD>
            y = testCase.y; X = testCase.X;
            testCase.startMeasuring();
            out = FSR(y, X, 'msg', false, 'plots', 0, 'nsamp', 500);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.beta);
        end

        function testFSRmdr(testCase, DataSize) %#ok<INUSD>
            y = testCase.y; X = testCase.X; bsb = testCase.bsb;
            testCase.startMeasuring();
            FSRmdr(y, X, bsb, 'msg', 0, 'plots', 0);
            testCase.stopMeasuring();
        end

        function testSreg(testCase, DataSize) %#ok<INUSD>
            y = testCase.y; X = testCase.X;
            testCase.startMeasuring();
            out = Sreg(y, X, 'msg', 0, 'plots', 0, 'nsamp', 500);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.beta);
        end

        function testMMreg(testCase, DataSize) %#ok<INUSD>
            y = testCase.y; X = testCase.X;
            testCase.startMeasuring();
            out = MMreg(y, X, 'plots', 0, 'Smsg', 0, 'nocheck', true);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.beta);
        end

        function testLXS_LMS(testCase, DataSize) %#ok<INUSD>
            y = testCase.y; X = testCase.X;
            testCase.startMeasuring();
            out = LXS(y, X, 'lms', 1, 'msg', false, 'plots', 0, 'nsamp', 1000);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.beta);
        end

        function testLXS_LTS(testCase, DataSize) %#ok<INUSD>
            y = testCase.y; X = testCase.X;
            testCase.startMeasuring();
            out = LXS(y, X, 'lms', 2, 'msg', false, 'plots', 0, 'nsamp', 1000);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.beta);
        end

        %% --- Multivariate ---

        function testFSM(testCase, DataSize) %#ok<INUSD>
            Y = testCase.Y;
            testCase.startMeasuring();
            out = FSM(Y, 'msg', false, 'plots', 0);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.md);
        end

        function testFSMmmd(testCase, DataSize) %#ok<INUSD>
            Y = testCase.Y; bsb = testCase.bsb;
            testCase.startMeasuring();
            FSMmmd(Y, bsb, 'msg', 0, 'plots', 0);
            testCase.stopMeasuring();
        end

        function testSmult(testCase, DataSize) %#ok<INUSD>
            Y = testCase.Y;
            testCase.startMeasuring();
            out = Smult(Y, 'msg', 0, 'plots', 0, 'nsamp', 200);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.loc);
        end

        function testMMmult(testCase, DataSize) %#ok<INUSD>
            Y = testCase.Y;
            testCase.startMeasuring();
            out = MMmult(Y, 'plots', 0, 'nocheck', true);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.loc);
        end

        function testMcd(testCase, DataSize) %#ok<INUSD>
            Y = testCase.Y;
            testCase.startMeasuring();
            out = mcd(Y, 'msg', 0, 'plots', 0);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.loc);
        end

        function testMahalFS(testCase, DataSize) %#ok<INUSD>
            Y = testCase.Y;
            mu = testCase.mu;
            Sigma = testCase.Sigma;
            testCase.startMeasuring();
            mahalFS(Y, mu, Sigma);
            testCase.stopMeasuring();
        end

        %% --- Clustering ---

        function testTclust(testCase, DataSize)
            Y = testCase.Y;
            k = DataSize.k;
            testCase.startMeasuring();
            out = tclust(Y, k, 0.05, 12, 'msg', 0, 'plots', 0, 'nsamp', 200);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.idx);
        end

        function testTkmeans(testCase, DataSize)
            Y = testCase.Y;
            k = DataSize.k;
            testCase.startMeasuring();
            out = tkmeans(Y, k, 0.05, 'msg', 0, 'plots', 0, 'nsamp', 200);
            testCase.stopMeasuring();
            testCase.verifyNotEmpty(out.idx);
        end
    end
end
