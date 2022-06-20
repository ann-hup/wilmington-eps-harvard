function G=JuliaMRSTmapping(JuliaFile)
    %{
        This function loads the grid information from Julia and creates the
        grid structure object required by MRST
        
        Input: 
            JuliaFile: julia file containing the grid information
    
     %}

    %Loading Julia input
    %load('/Users/josimar/Juanes Group Dropbox/Josimar Alves da Silva/Work/Projects/Induced_Seismicity/ExxonMobil/Modeling/BenchmarkProblem/ReservoirSimulation/validation/3D/fromJulia/CellStruct.mat')
    load(JuliaFile)
    
    %Creat dummy grid object to be changed
    G = removeCells( cartGrid([3,2]), 2);
    
    %Updating the grid structure
    G.cells.num = double(max(neighbors(:)));
    G.cells.facePos = double(facePos);
    G.cells.faces = double(faces(:,2));
    G.cells.indexMap = double([1:G.cells.num]');
    G.cells.cellPoro=double(cellPoro);
    G.cells.cellPerm=double(cellPerm);

    G.faces.num = double(size(faceArea,1));  %Number of faces in the model (including exterior)
    G.faces.nodePos = double(nodePos);
    G.faces.neighbors = double(neighbors);  
    G.faces.nodes = double(nodes(:,2));
    G.faces.tag=double(zeros(G.faces.num,1));
    
    G.nodes.num = double(size(coord,1));
    G.nodes.coords = double(coord);

    G.griddim = double(size(coord,2));
    G.type={'tensorGrid'};
    G.cartDims=[];

    G=computeGeometry(G);
    
    
    G.cells.centroids = double(cellCentroid);
    G.cells.volumes = double(cellVolume);
    
    G.faces.areas = double(faceArea);
    G.faces.normals= double(faceNormal);
    G.faces.centroids= double(faceCentroid);
    
    G.faces.TransJS=double(T);   %% This the transmissibility computed in Julia. T1, T2 and T=T1*T2/(T1+T2)
    G.cells.materialId = materialId;
end
