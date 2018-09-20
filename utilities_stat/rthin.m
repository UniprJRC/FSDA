function [Y , retain]= rthin(X, P)
%rthin applies independent random thinning to a point pattern.
%
%<a href="matlab: docsearchFS('rthin')">Link to the help page for this function</a>
%
% This function was ported to matlab from the R spatstat package, developed
% by Adrian Baddeley (Adrian.Baddeley@curtin.edu.au), Rolf Turner
% (r.turner@auckland.ac.nz) and Ege Rubak (rubak@math.aau.dk) for the
% statistical analysis of spatial point patterns. The algorithm for random
% thinning was changed in spatstat version 1.42-3. Our matlab porting is
% based on a earlier version. See the rthin documentation in spatstat for
% more details.
%
% In a random thinning operation, each point of X is randomly either
% deleted or retained (i.e. not deleted). The result is a point pattern,
% consisting of those points of X that were retained. Independent random
% thinning means that the retention/deletion of each point is independent
% of other points.
%
%  Required input arguments:
%
%       X   : Vector with the data to be thinned. Data can represent a
%             point pattern.
%
%       P   : Vector giving the retention probabilities, i.e. the probability
%             that each point in X will be retained. It can be:
%             -  a single number, so that each point will be retained with
%                the same probability P;
%             -  a vector of numbers, so that the ith point of X will be
%                retained with probability P(i);
%             -  a function P(x,y), so that a point at a location (x,y)
%                will be retained with probability P(x,y);
%             -  a pixel image, containing values of the retention
%                probability for all locations in a region encompassing the
%                point pattern.
%
%             If P is a function, it should be vectorised, that is, it
%             should accept vector arguments x,y and should yield a numeric
%             vector of the same length. The function may have extra
%             arguments which are passed through the argument.
%
%
%  Optional input arguments:
%
%
%  Output:
%
%    Y    :  the retained data units. Vector. In practice, Y = X(retain,:).
%
%  retain :  the indices of the retained points in the original data X.
%            Vector. The ith point of X is retained with probability P(i).
%
%
%  Optional Output:
%
%
%
%  See also ksdensity
%
%
% References:
%
% Bowman, A.W. and Azzalini, A. (1997), "Applied Smoothing
% Techniques for Data Analysis", Oxford University Press.
%
%
%
% Acknowledgements: 
%
% This function was ported to matlab from the R spatstat
% package, developed by Adrian Baddeley (Adrian.Baddeley@curtin.edu.au),
% Rolf Turner (r.turner@auckland.ac.nz) and Ege Rubak (rubak@math.aau.dk)
% for the statistical analysis of spatial point patterns. The algorithm for
% random thinning was changed in spatstat version 1.42-3. Our matlab
% porting is based on a previous version. See the rthin documentation in
% spatstat for more details.
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('rthin')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:
%
%{
    %% Random thinning on a mixture of normal distribution.
    % Data
    clear all; close all;
    data=[randn(500,2);randn(500,1)+3.5, randn(500,1);];
    x = data(:,1);
    y = data(:,2);

    % Data density
    [density,xout,bandwidth]   = kdebiv(data,'pdfmethod','fsda');
    xx = xout(:,1);
    yy = xout(:,2);
    zz = density;

    % plot of data and density
    figure;
    [xq,yq] = meshgrid(xx,yy);
    density = griddata(xx,yy,density,xq,yq);
    contour3(xq,yq,density,50), hold on
    plot(x,y,'r.','MarkerSize',5)
    title(['Original data (' num2str(numel(y)) ' units) with density contour'],'FontSize',16);

    %Interpolate the density and apply thinning using retention probabilities (1 - pdfe/max(pdfe))
    F = TriScatteredInterp(xx(:),yy(:),zz(:));
    pdfe = F(x,y);
    pretain = 1 - pdfe/max(pdfe);
    [Xt , Xti]= rthin([x y],pretain);

    % rthin retention probabilities
    [psorted ii] = sort(pretain);
    figure;
    plot(x,y,'r.','MarkerSize',5);
    hold on;
    plot(x(ii(1:100)),y(ii(1:100)),'bx','MarkerSize',5);
    title('The 100 units with smaller retention probabilities','FontSize',16);

    % now estimate the density on the retained units
    %[tdensity,txout,tbandwidth] = ksdensity(Xt);
    [tdensity,txout,tbandwidth]  = kdebiv(Xt,'pdfmethod','fsda');
    txx = txout(:,1);
    tyy = txout(:,2);
    tzz = tdensity;

    % and plot the retained units with their density superimposed
    figure;
    [txq,tyq] = meshgrid(txx,tyy);
    tdensity = griddata(txx,tyy,tdensity,txq,tyq);
    contour3(txq,tyq,tdensity,50), hold on
    plot(x(Xti),y(Xti),'b.','MarkerSize',5);
    title(['Retained data (' num2str(numel(y(Xti))) ' units) with new density contour'],'FontSize',16);

    cascade;
%}

%{
    %% Random thinning on the fishery dataset.
    % load data and add some jittering, because duplicated units are not treated
    clear all; close all;
    load('fishery.txt');
    fishery = fishery + 10^(-8) * abs(randn(677,2));
    x = fishery(:,1);
    y = fishery(:,2);

    % Data density
    [density,xout,bandwidth]   = kdebiv(fishery,'pdfmethod','fsda');
    xx = xout(:,1);
    yy = xout(:,2);
    zz = density;

    % plot of data and density
    figure;
    [xq,yq] = meshgrid(xx,yy);
    density = griddata(xx,yy,density,xq,yq);
    contour3(xq,yq,density,50), hold on
    plot(x,y,'r.','MarkerSize',8)
    xlim([0 300]); ylim([0 2000]);
    set(gca,'CameraPosition',[-216 -12425 0.0135]);
    title({['Zoom on fishery data (' num2str(numel(y)) ' units) with density contour'] , 'Probability mass concentrated close to the origin'},'FontSize',16);

    %Interpolate the density and apply thinning using retention
    %probabilities equal to 1 - pdfe/max(pdfe)
    F = TriScatteredInterp(xx(:),yy(:),zz(:));
    pdfe = F(x,y);
    pretain = 1 - pdfe/max(pdfe);
    [Xt , Xti]= rthin([x y],pretain);

    % now estimate the density on the retained units
    [tdensity,txout,tbandwidth]  = kdebiv(Xt,'pdfmethod','fsda');
    txx = txout(:,1);
    tyy = txout(:,2);
    tzz = tdensity;

    % and plot the retained units with their density superimposed
    figure;
    [txq,tyq] = meshgrid(txx,tyy);
    tdensity = griddata(txx,tyy,tdensity,txq,tyq);
    contour3(txq,tyq,tdensity,50), hold on
    plot(x(Xti),y(Xti),'b.','MarkerSize',8);
    xlim([0 300]); ylim([0 2000]);
    set(gca,'CameraPosition',[-216 -12425 0.0002558 ]);
    title({['Zoom on retained on the fishery data (' num2str(numel(y(Xti))) ' units) with density contour'] , 'Probabiity mass is smoother'},'FontSize',16);

    cascade;
 
%}


%% Beginning of code

n = length(X);

% if the retention probabilities are not provided by the user, retain
% with uniform probability 1/n
if nargin < 2
    P = 1/n;
end

if isnumeric(P)
    % vector of retention probabilities
    pX = P;
    if(length(pX) ~= n)
        if(length(pX) == 1)
            pX = repmat(pX, n, 1);
        else
            disp('Length of vector P does not match number of points of X');
            return;
        end
        if(any(isnan(pX)))
            disp('P contains NANs');
            return;
        end
    end
elseif isfunction(P)
    % function - evaluate it at points of X
    evalstr = [P '(x,y);'];
    pX = eval(evalstr);
    if(length(pX) ~= n)
        disp('Function P returned a vector of incorrect length');
        return;
    end
    if(~isnumeric(PX))
        disp('Function P returned non-numeric values');
    end
    if(any(isnan(pX)))
        disp('Function P returned some NA values');
        return;
    end
    %     rangemin = min(pX);
    %     rangemax = max(pX);
    %     prange = [min max];
end

if(min(pX) < 0)
    disp('some probabilities are negative');
    return;
end
if(max(pX) > 1)
    disp('some probabilities are greater than 1');
    return;
end

retain = rand(length(pX),1) < pX;

Y = X(retain,:);

%% this function is ported from the R function rthin from spatstat

% Acknowledgements: This function was ported to matlab from the R spatstat
% package, developed by Adrian Baddeley (Adrian.Baddeley@curtin.edu.au),
% Rolf Turner (r.turner@auckland.ac.nz) and Ege Rubak (rubak@math.aau.dk)
% for the statistical analysis of spatial point patterns. The algorithm for
% random thinning was changed in spatstat version 1.42-3. Our matlab
% porting is based on a previous version. See the rthin documentation in
% spatstat for more details.

%{
    % Skip_example
    % this is the version ported into matlab

    rthin <- function(X, P, ...) {
      verifyclass(X, "ppp")

      if(is.numeric(P)) {
        # vector of retention probabilities
        pX <- P
        if(length(pX) != X$n) {
          if(length(pX) == 1)
            pX <- rep(pX, X$n)
          else
            stop("Length of vector P does not match number of points of X")
        }
        if(any(is.na(pX)))
          stop("P contains NA's")
      } else if(is.function(P)) {
        # function - evaluate it at points of X
        pX <- P(X$x, X$y, ...)
        if(length(pX) != X$n)
          stop("Function P returned a vector of incorrect length")
        if(!is.numeric(pX))
          stop("Function P returned non-numeric values")
        if(any(is.na(pX)))
          stop("Function P returned some NA values")
        prange <- range(pX)
      } else if(is.im(P)) {
        # image - look it up
        if(!(P$type %in% c("integer", "real")))
          stop("Values of image P should be numeric")
        pX <- P[X, drop=FALSE]
        if(any(is.na(pX)))
          stop("some points of X lie outside the domain of image P")
      } else
      stop("Unrecognised format for P")

      if(min(pX) < 0) stop("some probabilities are negative")
      if(max(pX) > 1) stop("some probabilities are greater than 1")

      retain <- (runif(length(pX)) < pX)

      Y <- X[retain]

      # also handle offspring-to-parent map if present
      if(!is.null(parentid <- attr(X, "parentid")))
        attr(Y, "parentid") <- parentid[retain]

      return(Y)
    }
%}


%{
    % Skip_example
    % this is a newer version: to be  checked for consistency with older one

function (X, P, ..., nsim = 1, drop = TRUE)
    {
        verifyclass(X, "ppp")
        nX <- npoints(X)
        if (nX == 0) {
            if (nsim == 1 && drop)
                return(X)
            result <- rep(list(X), nsim)
            names(result) <- paste("Simulation", 1:nsim)
            return(as.solist(result))
        }
        if (is.numeric(P) && length(P) == 1 && spatstat.options("fastthin")) {
            result <- vector(mode = "list", length = nsim)
            for (isim in 1:nsim) {
                retain <- thinjump(nX, P)
                Y <- X[retain]
                if (!is.null(parentid <- attr(X, "parentid")))
                    attr(Y, "parentid") <- parentid[retain]
                result[[isim]] <- Y
            }
            if (nsim == 1 && drop)
                result <- result[[1]]
            return(result)
        }
        if (is.numeric(P)) {
            pX <- P
            if (length(pX) != nX) {
                if (length(pX) == 1)
                    pX <- rep.int(pX, nX)
                else stop("Length of vector P does not match number of points of X")
            }
            if (anyNA(pX))
                stop("P contains NA's")
        }
        else if (is.function(P)) {
            pX <- P(X$x, X$y, ...)
            if (length(pX) != nX)
                stop("Function P returned a vector of incorrect length")
            if (!is.numeric(pX))
                stop("Function P returned non-numeric values")
            if (anyNA(pX))
                stop("Function P returned some NA values")
        }
        else if (is.im(P)) {
            if (!(P$type %in% c("integer", "real")))
                stop("Values of image P should be numeric")
            pX <- P[X, drop = FALSE]
            if (anyNA(pX))
                stop("some points of X lie outside the domain of image P")
        }
        else stop("Unrecognised format for P")
        if (min(pX) < 0)
            stop("some probabilities are negative")
        if (max(pX) > 1)
            stop("some probabilities are greater than 1")
        if (nsim == 1) {
            retain <- (runif(length(pX)) < pX)
            Y <- X[retain]
            if (!is.null(parentid <- attr(X, "parentid")))
                attr(Y, "parentid") <- parentid[retain]
            return(if (drop) Y else solist(Y))
        }
        result <- vector(mode = "list", length = nsim)
        for (isim in 1:nsim) {
            retain <- (runif(length(pX)) < pX)
            Y <- X[retain]
            if (!is.null(parentid <- attr(X, "parentid")))
                attr(Y, "parentid") <- parentid[retain]
            result[[isim]] <- Y
        }
        names(result) <- paste("Simulation", 1:nsim)
        return(as.solist(result))
    }


%}

end

%FScategory:UTISTAT
