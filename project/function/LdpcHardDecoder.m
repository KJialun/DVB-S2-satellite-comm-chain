function [decoded_bits] = LdpcHardDecoder(bits_rx, H, numIt)

[m,n] = size(H); % row =m  collumn =n  
% the 1s of each raw represent the connection relationship between checknode to varible node

for i=1:m    % i iterate each check node
    % first we need to find out the index of each row non-zero elements (the variable nodes connectde to this checknode)
    c_to_v_nodes = find(H(i,:));
    for j=1:length(c_to_v_nodes)  %% j iterate each variable node of check node 
        c_nodes(i).msg = 0; % msg sent by c_node
        c_nodes(i).v_nodes(j).num = c_to_v_nodes(j); % index of each variable node connected to this checknode
        c_nodes(i).v_nodes(j).msg = 0; % msg received by c_node
    end
end
    
for i=1:n   % i iterate each variable node
    v_to_c_nodes = find(H(:,i));  % the index of check nodes  connected to this variable node
    for j=1:length(v_to_c_nodes)
        v_nodes(i).msg = 0; % msg sent by v_node
        v_nodes(i).c_nodes(j).num = v_to_c_nodes(j);
        v_nodes(i).c_nodes(j).msg = 0; % msg received by v_node
    end
end

graph = struct();
graph.c_nodes = c_nodes;
graph.v_nodes = v_nodes;


[M,N] = size(H);
it = 0;
decoded_bits = bits_rx;
% first step: 
% initialize the v_nodes
for i=1:N
    graph.v_nodes(i).msg = bits_rx(i);
end
% then communicate the msg to c_nodes
while ((it < numIt) && (length(find(mod(decoded_bits*H',2))) ~= 0))
    % communicate the msg to c_nodes
    for i=1:N
        c_nodes = graph.v_nodes(i).c_nodes; % c_nodes connected to a given v_node
        for j=1:length(c_nodes)
            idx = c_nodes(j).num; % get the global c_node index
            nums = [graph.c_nodes(idx).v_nodes(:).num];  % v_nodes connected the c_node
            for k=1:length(nums)
                % to get relative index of v_node with respect to c_node
                if nums(k) == i
                    graph.c_nodes(idx).v_nodes(k).msg = graph.v_nodes(i).msg;
                end
            end
        end
    end
    % second step:
    % calculate the response of the c_nodes using parity check equation
    for i=1:M
        v_nodes = graph.c_nodes(i).v_nodes; % v_nodes connected to a given c_node
        sumMsg = sum([v_nodes(:).msg]);
        for j=1:length(v_nodes)
            idx = v_nodes(j).num;
            nums = [graph.v_nodes(idx).c_nodes(:).num]; % c_nodes connected the v_node
            parityCheckMsg = mod(sumMsg - v_nodes(j).msg, 2);
            for k=1:length(nums)
                % to get relative index of v_node with respect to c_node
                if nums(k) == i
                    graph.v_nodes(idx).c_nodes(k).msg = parityCheckMsg;
                end
            end
        end
    end

    % third step:
    % v_nodes take a decision based on the response from c_node and original
    % msg using a majoroty vote
    for i=1:N
        c_nodes = graph.v_nodes(i).c_nodes(:);
        msgVector = [c_nodes.msg graph.v_nodes(i).msg];
        numVote = length(msgVector);
        voteBit1 = length(find(msgVector));
        percBit1 = voteBit1/numVote;
        if percBit1 > 0.5
            graph.v_nodes(i).msg = 1;
        else
            graph.v_nodes(i).msg = 0;
        end
    end
    
    decoded_bits = [graph.v_nodes(:).msg];
    it = it + 1;
end

end