function Vectors = PIV_conversion(PIVParams, Vectors)

%
% PIV_conversion - convert the PIV-vectors to physical values
% (This function is called by the function PIV)
%

Vectors.u =  Vectors.dx * PIVParams.FPS * PIVParams.Scaling / PIVParams.CorrelationStep;
Vectors.v = -Vectors.dy * PIVParams.FPS * PIVParams.Scaling / PIVParams.CorrelationStep;