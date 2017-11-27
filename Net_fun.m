classdef Net_fun < handle
    % Writen by winter-my-dream@hotmail.com
    % a program for generating defferent network
    methods (Static = true)
        %% NW-smallworld
        function matrix = small_world_NW(N,m,p)
            % 基于Matlab 的NW小世界网络仿真
            % 经过矩阵化修改后，生成速度已经大大加快
            % example: N=50;m=3;p=0.1;
            matrix = zeros(N,N);
            for i=m+1:N-m
                matrix(i,i-m:i+m)=1;
            end
            for i=1:m
                matrix(i,1:i+m)=1;
            end
            for i=N- m+1:N
                matrix(i,i-m:N)=1;
            end
            for i=1:m
                matrix(i,N-m+i:N)=1;
                matrix(N- m+i:N,i)=1;
            end
            % Randomly add edge
            kk=(rand(N,N)<p);
            matrix = logical(matrix + kk);
            matrix = matrix -diag(diag(matrix));
        end
        %% WS-smallworld
        function matrix = small_world_WS(N,m,p)
            % cite: Watts D J, Strogatz S H. Collective dynamics of ‘small-world’
            % networks[J]. nature, 1998, 393(6684): 440-442.
            % The simulation of WS-smallworld network
            % the algorithm of WS-smallworld's generation has been improved in speed,
            % and tend to be easily understood
            % Example:
            % N = 100; %network size (number of nodes)
            % m = 6; %2*m is the average edges of each nodes
            % p = 0.1; %rewiring probability
            % matrix = small_world_WS_new(N,m,p);
            rng('default')
            rng('shuffle')
            matrix=zeros(N,N);
            % generate regular network
            for i=m+1:N-m
                matrix(i,i-m:i+m)=1;
            end
            for i=1:m
                matrix(i,1:i+m)=1;
            end
            for i=N-m+1:N
                matrix(i,i-m:N)=1;
            end
            for i=1:m
                matrix(i,N-m+i:N)=1;
                matrix(N-m+i:N,i)=1;
            end
            % rewiring the network
            for i = 1:N
                % then rewiring the edges with the probability of p
                [series1,series2] = range_sort(N,m,i);
                index0 = series1(rand(2*m,1)>1-p);
                if(~isempty(index0))
                    matrix(i,index0) = 0;
                    matrix(i,series2(randperm(length(series2),length(index0))))=1;
                end
            end
            matrix = matrix -diag(diag(matrix));
        end
        function [series1,series2] = range_sort(N,m,i)
            % select the index of nodes in row i for rewiring
            if(i-m>0 && i+m<=N)
                series1 = i-m:i+m;
                series2 = setdiff(1:N,series1);
            elseif(i-m<=0)
                series1 = [1:i+m,N-m+i:N];
                series2 = setdiff(1:N,series1);
            else
                series1 = [1:m-N+i,i-m:N];
                series2 = setdiff(1:N,series1);
            end
            % Without considering the connection of diagonal elements
            series1(series1==i) = [];
        end
        %% scale-free network
        function a=sfnet(N,p)
            % Generate a scale free network using the BA algorithm
            % N = number of nodes
            % p <-> number of links
            % By default, the network is symmetrical (a_ij=a_ji)
            % but columns and/or lines can be randomized
            % Example: a=SF_net(50,0.5)
            
            % Adapted by Didier Gonze
            % Created: 25/4/2014
            % Updated: 25/4/2014
            
            seed =[0 1 0 0 1;1 0 0 1 0;0 0 0 1 0;0 1 1 0 0;1 0 0 0 0];
            m=round(p*N/10); % number of links
            
            a = SFNG(N,m,seed);
            
            % Randomization:
            
            %a=a(:,randperm(size(a,1)));  % randomize lines
            %a=a(:,randperm(size(a,2)));  % randomize columns
            
            a=double(a);
        end
        % ============================================================
        % Function SFNG
        % ============================================================
        % Function taken from:
        % http://www.mathworks.com/matlabcentral/fileexchange/11947-b-a-
        % scale-free-network-generation-and-visualization/content/SFNG.m
        function Net = SFNG(Nodes, mlinks, seed)
            seed = full(seed);
            pos = length(seed);
            
            %if (Nodes < pos) || (mlinks > pos) || (pos < 1) || (size(size(seed)) ~= 2) || (mlinks < 1) || (seed ~= seed') || (sum(diag(seed)) ~= 0)
            %    error('invalid parameter value(s)');
            %end
            
            %if mlinks > 5 || Nodes > 15000 || pos > 15000
            %    warning('Abnormally large value(s) may cause long processing time');
            %end
            rand('state',sum(100*clock));
            Net = zeros(Nodes, Nodes, 'single');
            Net(1:pos,1:pos) = seed;
            sumlinks = sum(sum(Net));
            
            while pos < Nodes
                pos = pos + 1;
                linkage = 0;
                while linkage ~= mlinks
                    rnode = ceil(rand * pos);
                    deg = sum(Net(:,rnode)) * 2;
                    rlink = rand * 1;
                    if rlink < deg / sumlinks && Net(pos,rnode) ~= 1 && Net(rnode,pos) ~= 1
                        Net(pos,rnode) = 1;
                        Net(rnode,pos) = 1;
                        linkage = linkage + 1;
                        sumlinks = sumlinks + 2;
                    end
                end
            end
            clear Nodes deg linkage pos rlink rnode sumlinks mlinks
        end
        %% Example of scale-free network
        function Scale_free_examplae
            tic
            m0 = 10;
            m = 5;
            NodesNum = 1000;
            A = sparse(NodesNum, NodesNum);
            A(1:m0,1:m0) = round(rand(m0));
            A = tril(A);
            A = A+A';
            A = A - diag(diag(A));
            for i=m0+1:NodesNum
                Degree = sum(A(1:i-1,1:i-1));
                for j=2:i-1
                    Degree(j) = Degree(j) + Degree(j-1);
                end
                LinksNum = 0;
                while LinksNum<m
                    link = fix(rand()*Degree(i-1)+1);
                    for j=1:i-2
                        if link<=Degree(1) && A(i,1)==0
                            A(i,1) = 1;
                            A(1,i) = 1;
                            LinksNum = LinksNum+1;
                        elseif link>Degree(j) && link<=Degree(j+1) && A(i,j+1)==0
                            A(i,j+1) = 1;       A(j+1,i) = 1;
                            LinksNum = LinksNum+1;
                        end
                    end
                end
            end
            Degree = sum(A);
            list = unique(Degree);
            num = zeros(1,length(list));
            for i=1:length(list)
                num(i) = length(find(list(i)==Degree));
            end
            toc
            loglog(list,num ./ sum(num),'.','markersize',20)
            xlabel('k'),ylabel('P(k)')
        end
        %% visualize network
        function plot_network(N,m,p,type)
            clear,clc,close all
            % example:
            %             N=10;
            %             m=2;
            %             p=0.1;
            switch type
                case 'ws-smallworld'
                    A= small_world_WS_new(N,m,p);
                case 'nw-smallworld'
                    A = small_world_NW(N, m, p);
            end
            t=linspace(0,2*pi,N+1);
            x=sin(t);
            y=cos(t);
            figure
            set(gcf,'color','w')
            plot(x,y,'o','markerfacecolor','k'),hold on
            for i=1:N
                for j=1:N
                    if (A(i,j)==1)
                        fp1=plot([x(i),x(j)],[y(i),y(j)],'r-'); hold on
                        set(fp1,'linesmoothing','on')
                    end
                end
            end
            axis([-1.05,1.05,-1.05,1.05])
            axis square
            axis off
            sum(sum(A))
        end
    end
end
