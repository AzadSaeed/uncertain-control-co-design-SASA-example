function X = BrysonDenham_Y(t,l)

    X = zeros(length(t),2);
    
    %----------------------------------------------------------------------
    if (l > 0) && (l <= 1/6)
        for i = 1:length(t)
            if (t(i) <= 3*l)
                v = (1 - t(i)/(3*l)).^2;
                x = l*(1 - (1 - t(i)/(3*(l))).^3);
            elseif (t(i) >= 1-3*l)
                v = -(1 - (1-t(i))/(3*l)).^2;
                x = l*(1 - (1 - (1-t(i))/(3*(l))).^3);
            else
                v = 0;
                x = l;
            end
            % create the state matrix
            X(i,:) = [x,v];
        end
    %----------------------------------------------------------------------
    elseif (l >= 1/6) && (l <= 1/4)
        for i = 1:length(t)
            if (t(i) <= 1/2)
                v = 1 - 8*(1-3*l)*t(i) + 12*(1-4*l)*(t(i)).^2;
                x = t(i) - 4*(1-3*l)*(t(i)).^2 + 4*(1-4*l)*(t(i)).^3;
            else
                v = -1 + 8*(1-3*l)*(1-t(i)) - 12*(1-4*l)*(1-t(i)).^2;
                x = 1 - t(i) - 4*(1-3*l)*(1-t(i)).^2 + 4*(1-4*l)*(1-t(i)).^3;
            end
            % create the state matrix
            X(i,:) = [x,v];
        end
    %----------------------------------------------------------------------
    elseif (l >= 1/4)
        for i = 1:length(t)
            v = 1-2*t(i);
            x = t(i)*(1-t(i));
            % create the state matrix
            X(i,:) = [x,v];
        end
    %----------------------------------------------------------------------
    else
        X = nan(length(t),2);
    end

end