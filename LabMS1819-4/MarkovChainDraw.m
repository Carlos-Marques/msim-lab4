

nodes=20;   
%length of a side of square in which sensors are located
sidelength=100;

% generating node id and x,y positions
nodePos=[ (1:nodes)' (rand(nodes,2)*sidelength) ];

% loading the saved positions
load MarkovChain
nodes=20;
figure
%plotting the nodes along with their id
plot(nodePos(:,2),nodePos(:,3),'ob');
text(nodePos(:,2),nodePos(:,3),num2str(nodePos(:,1)));
hold on

% list of edges
node_list=[15 5;5 11;11 6;15 6;1 20;1 7; 7 20;7 16; 16 18;18 14 ;14 20; 20 1;20 7; 1 6;7 19;13 19;2 13;2 4;13 4; 19 3; 3 12 ;12 10;12 8; 8 9;10 9;10 17; 9 17; 4 19];

%creating adjacency matrix
A=zeros(nodes,nodes);
[mtrue,~]=size(node_list);
for create=1:mtrue
A(node_list(create,1),node_list(create,2))=1;
A(node_list(create,2),node_list(create,1))=1;
end

%plotting an edge between two nodes if they are connected by an edge
for k=1:nodes
  for i=1:nodes
       if A(k,i)~=0
              
      line=[nodePos(k,2:3);nodePos(i,2:3)];
     plot(line(:,1),line(:,2))
      end
  end
end

%plotting the source position
plot(sourcePos(:,1),sourcePos(:,2),'rh');