route_in = '/Users/Endao/Desktop/Ultrasound/US_Code/PIV_Parameters';
route_out = '/Volumes/Ultrasound/Ultrasound/20150427';

for vp = 50:50:50
    for i = 1:1:1
        
        name = ['cs_' num2str(vp) 'mms_' num2str(i) '_recon'];
        
        PIVParams = PIV_readPIVParameters([route_in '/cs' num2str(vp) 'mms.txt']);
        PIVParams.ImageFileName = 'bla';
        PIVParams.ROITopEdge = 91;
        PIV_writePIVParameters(PIVParams, [route_out '/' name '/Parameters.txt']);
        
    end
end