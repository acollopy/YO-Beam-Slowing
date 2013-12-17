function [XV,LostToDelta] = LossIntoDeltaAndv2(XV,vStart,StateStart,angle,v2BeamSize,LostToDelta)  
    NumberOfMolecules=size(XV,2);
    CellToMotDist = .44;
    
    %state: -1 for delta, 0 for v=2, 1 for v=0/1
      
    rsn = rand(1,NumberOfMolecules);
    LossFractionv2 = 1*( 1-exp(-(vStart-XV(6,:))/(6e-3*3000)));%Loss into v = 2;  3000 photons before going into v = 2, 6mm/s per recoil;  
    LossFractionDelta = 1*(1-exp(-(vStart-XV(6,:))/(6e-3*3000)));%Loss into Delta;  3000 photons before going into delta;  

    InV2 = zeros(1,NumberOfMolecules) + ((XV(1,:)-tan(angle)*(XV(3,:)-.5*CellToMotDist)).^2+(XV(2,:).^2) <= v2BeamSize^2);
    StateNew = StateStart | InV2; %if in v=2 beam, move back to state 1

    XV(9,:) = (ceil(rsn-LossFractionv2-LossFractionDelta) - ceil(LossFractionDelta-rsn)).*StateNew(1,:);
    a = size(XV,2);
    XV = XV(:,XV(9,:)>-1);%Remove Delta State
    LostToDelta = LostToDelta +a-size(XV,2);
end
