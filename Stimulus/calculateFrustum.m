function calculateFrustum(coordinateMuilty,deltaDegree)
global TRIALINFO;
global FRUSTUM;
global SCREEN;

FRUSTUM.clipNear = 50*coordinateMuilty; % cm
FRUSTUM.clipFar = 150*coordinateMuilty; % cm
FRUSTUM.top = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.heightCM / 2.0);
FRUSTUM.bottom = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.heightCM / 2.0);

if nargin ==1
    % left eye
    FRUSTUM.sinisterRight = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.widthCM / 2.0 + TRIALINFO.deviation / 2.0);
    FRUSTUM.sinisterLeft = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.widthCM / 2.0 + TRIALINFO.deviation / 2.0);
    
    % right eye
    FRUSTUM.dexterRight = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.widthCM / 2.0 - TRIALINFO.deviation / 2.0);
    FRUSTUM.dexterLeft = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.widthCM / 2.0 - TRIALINFO.deviation / 2.0);
    
elseif nargin == 2
    delta = FRUSTUM.clipNear * sind(deltaDegree);
    % left eye
    FRUSTUM.sinisterRight = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.widthCM / 2.0 + TRIALINFO.deviation / 2.0)+delta;
    FRUSTUM.sinisterLeft = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.widthCM / 2.0 + TRIALINFO.deviation / 2.0)+delta;
    
    % right eye
    FRUSTUM.dexterRight = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.widthCM / 2.0 - TRIALINFO.deviation / 2.0)+delta;
    FRUSTUM.dexterLeft = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.widthCM / 2.0 - TRIALINFO.deviation / 2.0)+delta;
end