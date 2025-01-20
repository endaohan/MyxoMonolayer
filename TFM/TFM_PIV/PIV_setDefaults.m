function PIVParams = PIV_setDefaults()

% Tasks
PIVParams.DoPIVAnalysis = true;
PIVParams.ConvertData = true;
PIVParams.CalcScalars = true;
PIVParams.SaveData = true;
PIVParams.DisplayEveryFrame = true;
PIVParams.SaveImages = true;

% File locations and names
PIVParams.ImageFileName = '';
PIVParams.FileExtension = '';
PIVParams.MultipageTiff = true;
PIVParams.ImageBitDepth = 8;
PIVParams.Directory = '';
PIVParams.ImageOutputDir = '';
PIVParams.ImageOutputFormat = 'png';
PIVParams.DataDir = '';
PIVParams.FirstImage = 1;
PIVParams.LastImage = 1;
PIVParams.SkipImages = 0;

% Image pre-processing
PIVParams.UseMinMax = false;
PIVParams.MinMaxKernel = 5;
PIVParams.UseUnsharpMask = false;
PIVParams.UnsharpMaskAmount = 1;
PIVParams.UseGaussianBlur = false;
PIVParams.GaussianBlurRadius = 0.3;

% Grid and ROI
PIVParams.RotationAngle = 0;
PIVParams.GridSize = 32;
PIVParams.Overlap = 0;
PIVParams.UseROI = false;
PIVParams.CropImageToROI = false;
PIVParams.ROILeftEdge = 0;
PIVParams.ROIRightEdge = 0;
PIVParams.ROITopEdge = 0;
PIVParams.ROIBottomEdge = 0;
PIVParams.UseMovingROI = false;
PIVParams.MovingROILeftStart = 0;
PIVParams.MovingROIRightStart = 0;
PIVParams.MovingROIHorizontalVelocity = 0;
PIVParams.MovingROITopStart = 0;
PIVParams.MovingROIBottomStart = 0;
PIVParams.MovingROIVerticalVelocity = 0;

% Correlation
PIVParams.Method = 'XCORR';
PIVParams.NumPasses = 1;
PIVParams.UseGridRefinement = false;
PIVParams.InitialGridSize = 64;
PIVParams.CorrelationStep = 1;

% Validation and filtering
PIVParams.numValidationPasses = 5;
PIVParams.UseMaxDisplacement = true;
PIVParams.MaxDisplacement = 10;
PIVParams.UseMinCorrelationCoefficient = false;
PIVParams.MinCorrelationCoefficient = 0.5;
PIVParams.UseNormalizedMedianTest = true;
PIVParams.MedianThreshold = 4;
PIVParams.UseSecondPeak = false;
PIVParams.UseInterpolation = false;
PIVParams.UseDataSmoothing = false;
PIVParams.SmoothingFilterSize = 3;
PIVParams.SmoothingFilterType = 'gaussian';

% Conversion
PIVParams.FPS = 1000;
PIVParams.Scaling = 1.5e-5;
PIVParams.TimeUnit = 's';
PIVParams.LengthUnit = 'm';
PIVParams.VelocityUnit = 'm/s';

% Visualization
PIVParams.ShowVectors = true;
PIVParams.ColorCodeVectors = false;
PIVParams.VectorColorValidation = false;
PIVParams.VectorColormap = 'jet';
PIVParams.VectorColor = 'black';
PIVParams.VectorScale = 10;
PIVParams.ShowImage = false;
PIVParams.ShowGrid = false;
PIVParams.ShowContour = true;
PIVParams.ContourScalar = 'magnitude';
PIVParams.ContourColormap = 'jet';
PIVParams.ContourSteps = 64;
PIVParams.ContourLines = false;
PIVParams.ScalarMinValue = 0;
PIVParams.ScalarMaxValue = 10;
PIVParams.ScalarAutoScale = false;
PIVParams.ScalarPrefactor = 1;
PIVParams.DisplayScalarInfo = true;
PIVParams.OutputImageWidth = 4;
PIVParams.OutputImageHeight = 3;
PIVParams.FontName = 'Helvetica';
PIVParams.FontSize = 10;
PIVParams.ScalarLabel = '';

% For internal use, not written to parameters file
PIVParams.CurrentParametersFile = '';
PIVParams.CurrentDir = cd;