%% Computes Weighted Contact Map, Pairwise Electrostatic Energy between Charged Atoms and Defined Block Approximation based on the Protein Size

%% Required Input Files - Modify them in the Input Parameters Section
% a)PDBfile - Make sure the pdb file is available in the same folder as
% this script. Also ensure there are no missing atoms/residues in which case the
% outputs will be erroneous
% b)Stride output file - Please use the same PDB file as above to obtain the results from STRIDE server
%                        URL - http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py
% c)pH - takes values of 5 or 7. At pH 5, His will be additionally charged apart from Lys, Arg, Asp and Glu
% d)Cutoff to identify VdW contacts (in Ang. units)

clear; clc; tic;

%% Input Parameters
pdblist = char('CI2'); % Input PDB filename
stridefile = char('struct.txt'); % File containing raw output from the STRIDE server
pH = 7; % pH value
srcutoff=5.0; % Cutoff distance for VdW interaction
BlockSize = 1; % can be manually set to a specific value. A heuristic for blocksize is also provided in lines 132-139


cdata=cat(2,deblank(pdblist(1,:)),'.pdb');
pdbr=fopen(cdata,'rt'); % Raw PDB data
pdbc=fopen('buffer1.txt','wt'); % Curated PDB data

%% Remove Hydrogens and hetero atoms including nucleotides
% gets the raw PDB file and generates the buffer1.txt required for
% calculations
line=fgetl(pdbr);
while feof(pdbr)~=1
    if(length(line)>4)
        if(strcmp(line(1:4),'ATOM')==1 && strcmp(line(13),'H')==0 && strcmp(line(13),'Q')==0 && strcmp(line(14),'H')==0 && strcmp(line(18),' ')==0)
            fprintf(pdbc,'%s\n',line);
        end
        if(strcmp(line(1:3),'TER')==1) % Terminates at the first instance - Please modify the input file accordingly if needed
            break;
        end
    end
    line=fgetl(pdbr);
end
fclose(pdbr);
fclose(pdbc);

%% Get atmoic details with their corresponding residue number
% atomn: Atom name, resno: Residue number, x,y,z: coordinates

pdbc=fopen('buffer1.txt','rt');
i=1; j=1; prevres=0; nnegres=0; nposres=0; nhphilres=0; iin=1; iim=1;
aares={'GLY','ALA','VAL','LEU','ILE','MET','PHE','TYR','TRP','SER','ASP','ASN','THR','GLU','GLN','HIS','LYS','ARG','PRO','CYS'};
aacode={'G','A','V','L','I','M','F','Y','W','S','D','N','T','E','Q','H','K','R','P','C'};
charres={'HIS','LYS','ARG','GLU','ASP'};
posres=[16 17 18]; % [HIS LYS ARG]
negres=[11 14]; % [ASP GLU]
hphilres=[8 10 12 13 15 20]; % [TYR SER ASN THR GLN CYS]
atomc={'NE ','NH1','NH2','NZ ','OD1','OD2','OE1','OE2','ND1','NE2'};
atombb={'C','N','CA'};
charmag7=[0.33 0.33 0.33 1 -0.5 -0.5 -0.5 -0.5 0 0]'; % pH 7
charmag5=[0.33 0.33 0.33 1 -0.5 -0.5 -0.5 -0.5 0.5 0.5]'; % pH 5

while feof(pdbc)~=1
    line=fgetl(pdbc);
    if(length(line)>4)
        atomn(i,:)=line(14:16);
        resno(i,1)=str2double(line(23:26));
        x(i,1)=str2double(line(32:38));
        y(i,1)=str2double(line(40:46));
        z(i,1)=str2double(line(48:54));
 
        % Assign charge magnitude at specified pH
        charmag(i,1)=0;
        if sum(strcmp(charres,line(18:20)))
            [~,x1] = find(strcmp(atomc,atomn(i,:)));
            if isempty(x1)==0
                if pH == 7
                    charmag(i,1)=charmag7(x1);
                elseif pH == 5
                    charmag(i,1)=charmag5(x1);
                end
            end
        end
        
        % Calculate Seq Composition
        if (resno(i)-prevres~=0)
            % Get sequence from PDB
            [~,x2]=find(strcmp(line(18:20),aares));
            protseq(1,j)=aacode(x2);
            prevres=resno(i);
            j=j+1;
        end
    end
    i=i+1;
end

%% Compute Short and Long-Range interaction map
% Short_Range: no. of atomic contacts and h-bonds b/w residues
% Long_Range: details for electrostatic interactions

ecutoff=1000; % Cutoff dist for electrostatic interaction
%-----Initializarion-----------------------------------------------------------------------
atomt=length(x(:,1));
nres=resno(atomt,1)-resno(1,1)+1;
srcont=zeros(nres);
%------------------------------------------------------------------------------------------
startres=resno(1,1);
resno(:,1)=resno(:,1)-startres+1; %Changing res. no. according to the PDB file
k=1;

for i=1:(atomt-1)
    for j=(i+1):atomt
        dist=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2+(z(i)-z(j))^2);
        % Calculate Short-Range Interactions
        if((resno(j)-resno(i))>=1 && (charmag(i)==0 || charmag(j)==0))
            dist=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2+(z(i)-z(j))^2);
            if(dist<=srcutoff)
                srcont(resno(i),resno(j))=srcont(resno(i),resno(j))+1;
            end
        end
        
        % Calculate Long-Range Interactions and associated Energy
        if ((resno(j)-resno(i))>=1 && (charmag(i)~=0 && charmag(j)~=0))
            diste=sqrt(((x(i)-x(j))^2)+((y(i)-y(j))^2)+((z(i)-z(j))^2));
            if(diste<=ecutoff)
                intene=332*4.184*charmag(i)*charmag(j)/(29*diste); %Elec. Contibution
                % Elec. Info in matrix form
                ElecMat(k,:)=[resno(i) resno(j) diste abs(resno(i)-resno(j)) intene];
                k=k+1;
            end
        end
    end
end

% Heuristic to choose the block size
% if nres<=300
%     BlockSize = ceil(nres/100);
% elseif nres<600
%     BlockSize = 5;
% else
%     BlockSize = 6;
% end

hh9=fopen(stridefile(1,:),'rt'); i=1; j=1;
lkk=fgetl(hh9);
while lkk > 0
    if length(lkk)> 25 && strcmp(lkk(1:3),'ASG')
        h(i,1) = 0;
        if (strcmp(lkk(25),'H') || strcmp(lkk(25),'E') || strcmp(lkk(25),'G'))
            h(i,1) = 1;
        end
        i=i+1;
    end
    lkk=fgetl(hh9);
end
fclose(hh9);

h1=[[0;h] [h;0]];
hbeg = find(h1(:,1)==0 & h1(:,2)==1);
hend = find(h1(:,1)==1 & h1(:,2)==0);
StructBlock = startres; hbeg = hbeg+startres-1; hend=hend+startres-1;
for i=1:length(hbeg)
    StructBlock = [StructBlock;hbeg(i);hend(i)];
end
StructBlock = [StructBlock;resno(end)+startres];
residual = rem(StructBlock(2:end)-StructBlock(1:end-1),BlockSize);
nBlocks = floor((StructBlock(2:end)-StructBlock(1:end-1))/BlockSize);

k=1; BlockID = 0; uresno = unique(resno);
BlockMat = zeros(length(uresno),2); flag=0;
for i=2:length(StructBlock)
    for j=1:nBlocks(i-1)
        BlockID = BlockID+1;
        BlockMat(k:k+BlockSize-1,:) = [uresno(k:k+BlockSize-1) BlockID*ones(BlockSize,1)];
        k=k+BlockSize; flag=1;
    end
    if residual(i-1) == 1
        if flag==0; BlockID = BlockID+1; end
        BlockMat(k,:) = [uresno(k) BlockID];
        k=k+1;
    elseif residual(i-1) > 1
        BlockID = BlockID+1;
        BlockMat(k:k+residual(i-1)-1,:) = [uresno(k:k+residual(i-1)-1) BlockID*ones(residual(i-1),1)];
        k=k+residual(i-1);
    end
end

contactMapB = zeros(BlockID);
for i=1:nres
    for j=i+1:nres
        contactMapB(BlockMat(i,2),BlockMat(j,2)) = contactMapB(BlockMat(i,2),BlockMat(j,2)) + srcont(i,j);
    end
end

contdistMapB = zeros(length(ElecMat(:,1)),5);
for i=1:length(ElecMat(:,1))
    contdistMapB(i,:) = [BlockMat(ElecMat(i,1),2) BlockMat(ElecMat(i,2),2) ElecMat(i,3:end)];
end

%% Printing Regular Comtact Map files
eval(['save contactmapmatElec',deblank(pdblist(1,:)),'.dat srcont -ascii']);
eval(['save contdistElec',deblank(pdblist(1,:)),'.dat ElecMat -ascii']);
eval(['save contdistElecB',deblank(pdblist(1,:)),'.dat contdistMapB -ascii']);
eval(['save contactmapmatElecB',deblank(pdblist(1,:)),'.dat contactMapB -ascii']);
eval(['save BlockDet',deblank(pdblist(1,:)),'.dat BlockMat -ascii']);
eval(['save BlockSize',deblank(pdblist(1,:)),'.dat BlockSize -ascii']);
