function [eps, hest] = SolveEpsilonCorrection(solver, sigmainit)



    iter = 1; %The number of current iteration
    maxIter = 30;
    erel = 1e-4;
    
    hfig = figure();

    res = solver.OptimizationFunction(0*sigmainit); %the residual of minimum conductivity

    cont = 1;%continue flag for the linesearch
    d = [0, 1]';%vector containing different relatve step lengths tried
    ii = 2;%the line search iteration counter
    r = 0; %the index of the right point (0 = we're still checking if we want to go even further right)
    l = 1; %the index of the left point
    c = 1; %the index of the center point
    while cont == 1
        sigest_new = d(ii)*sigmainit;
        %sigest_new(sigest_new < sigmamin) = sigmamin;
        res(ii) = solver.OptimizationFunction(sigest_new);%in each point we calculate the residual
        if r == 0 %we're still going right
            if res(ii) < res(ii-1) %continue right
                l = ii-1;
                d(ii+1) = 2*d(ii);
            else%stop going right
                r = ii;
                c = ii-1;
                d(ii+1) = 0.5*(d(r)+d(c));
            end
        else%the minimum is between our left and right points
            if (d(ii)<d(c) && res(ii) < res(c))
                r = c;
                c = ii;
                d(ii+1) = 0.5*(d(c)+d(r));
            elseif (d(ii)<d(c) && res(ii) > res(c))
                l = ii;
                d(ii+1) = 0.5*(d(c)+d(r));
            elseif (d(ii)>d(c) && res(ii) < res(c))
                l = c;
                c = ii;
                d(ii+1) = 0.5*(d(c)+d(r));
            elseif (d(ii)>d(c) && res(ii) > res(c))
                r = ii;
                if (c~=l)
                    d(ii+1) = 0.5*(d(l)+d(c));
                else
                    d(ii+1) = 0.5*(d(l)+d(ii));
                end
            else%we should not end up here
                warning('Something wrong with line search');
                break;
            end
        end

        set(0, 'CurrentFigure', hfig);plot(d(1:end-1), res, 'bo', d(ii), res(ii), 'r+'); %plot the progress of the linesearch
        hfig.Name = ['Iteration ', num2str(iter)];%mark down the iteration too on this window
        drawnow;
        
        if (r > 0)%see if we're happy with our results
            if (abs(res(l) - min(res)) < erel*res(ii) && abs(res(r) - min(res)) < erel*res(ii))
                cont = 0;
            end
        end
        ii = ii +  1;
        if (ii > maxIter)
            cont = 0;
        end

    end


    close(hfig);
    
    [~,minii] = min(res);
    hest = d(minii);
    elval = solver.SolveForwardVec(hest*sigmainit);
    eps = solver.Dmeas - elval;
    

end