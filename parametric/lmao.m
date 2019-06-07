clear all

% Missing file 1 (AR 20)
[x, Fs]=audioread('../audio/armst_37_orig.wav');
xt=x(:,1)';
error_var=9.4588e-08;

P=20;


mmselist=[];

for param=1:20
    N=size(xt,2);
    xt_lossy=xt;
    
    missing=param*10;
    anchors=linspace(200,N-200,18);
    for i=1:length(anchors)
        anchor=anchors(i);
        p=round(anchor-(missing/2-1));
        q=round(anchor+(missing/2));
        xt_lossy(p:q)=0;
    end

    x_f = xt_lossy;
    x_b = xt_lossy;
    x_predict1 = xt_lossy;

    checkarray=[1];
    counter=1;
    for i=1:N
        if counter==1
            if xt_lossy(i)==0
                checkarray=[checkarray i];
                counter=0;
            else
                continue
            end
        else
            if xt_lossy(i)==0
                continue
            else
                checkarray=[checkarray i-1];
                counter=1;
            end
        end
    end
    checkarray=[checkarray N];

    [b, i1] = unique(checkarray,'first');
    [b, i2] = unique(checkarray,'last');
    checkarray = b(i1==i2);

    o=checkarray(1);
    p=checkarray(2);
    q=checkarray(3);
    r=checkarray(4)-1;

    sizearray=size(checkarray,2);
    numZeros=0.5*(sizearray-2);
    checkarray=[checkarray 0 0];

    for i=1:numZeros 
        data=xt(o:p-1);
        N=size(data,2);
        G = fliplr(buffer(data(1:end-1), N-P, N-1-P, 'nodelay'));
        y=data(P+1:end)';
        % ML estimate for theta
        theta_ML = inv(transpose(G)*G)*transpose(G)*y;
        % Prior distribution parameters
        theta_prior = zeros(P,1);
        prior_var = eye(P);
        % Posterior parameters (inspecting values, closer to 0 due to the prior)
        [theta_MAP,phi,big_theta,post_var] = post_param(G,error_var,...
            theta_prior,prior_var,y);
        L=q-p+1;

        % Forward prediction mode (w/o adding noise)
        for packet = p:p+L-1
            x_f(packet) = x_f(packet-1:-1:packet-P)*theta_MAP;
        end

        % Backward prediction mode (w/o adding noise)
        for packet = p+L-1:-1:p
            x_b(packet) = x_b(packet+1:packet+P)*theta_MAP;
        end

    %     % Forward prediction mode (w adding noise)
    %     for packet = p:p+L-1
    %         x_f(packet) = x_f(packet-1:-1:packet-P)*theta_MAP + ...
    %             sqrt(error_var)*randn(1,1);
    %     end
    % 
    %     % Backward prediction mode (w adding noise)
    %     for packet = p+L-1:-1:p
    %         x_b(packet) = x_b(packet+1:packet+P)*theta_MAP + ...
    %             sqrt(error_var)*randn(1,1);
    %     end

        % Weighted sum
        for j = 1:L
            alpha = (L-j)/(L-1);
            x_predict1(p-1+j) = alpha*x_f(p-1+j) + (1-alpha)*x_b(p-1+j);
        end


        o=q+1;
        p=r+1;
        q=checkarray(3+i*2);
        r=checkarray(4+i*2)-1;
    end

    mmselist=[mmselist, immse(xt, x_predict1)];
end

x=1:20;
x=x*10;
plot(x, mmselist)
xlabel('Length of missing data sections')
ylabel('MSE in interpolation')
title('Mean-squared error in interpolation')