function tests = test_speed_reducer
% TEST_SPEED_REDUCER  Function based unit tests (matlab.unittest).
%
% Author: AWD Labs
% Student ID: 52104479
%
% Run from code/matlab:  results = runtests('test_speed_reducer');
tests = functiontests(localfunctions);
end

function xv = xval()
xv = [3.5; 7.0; 22; 7.4; 7.8; 3.5; 5.2];
end

function testStressesAt2400(tc)
d = speed_reducer_analysis(xval(), speed_reducer_params('validation2400')).display;
verifyEqual(tc, d.sigma_b_MPa, 323, 'AbsTol', 1);
verifyEqual(tc, d.sigma_c_MPa, 532, 'AbsTol', 1);
verifyEqual(tc, d.sigma_s_MPa(1), 512, 'AbsTol', 1);
verifyEqual(tc, d.sigma_s_MPa(2), 454, 'AbsTol', 1);
end

function testDeflectionsAt2400(tc)
d = speed_reducer_analysis(xval(), speed_reducer_params('validation2400')).display;
verifyEqual(tc, d.y_mm(1), 0.0179, 'AbsTol', 5e-5);
verifyEqual(tc, d.y_mm(2), 0.0043, 'AbsTol', 5e-5);
end

function testTorque1000DoesNotMatch(tc)
d = speed_reducer_analysis(xval(), speed_reducer_params('validation1000')).display;
verifyLessThan(tc, d.sigma_b_MPa, 0.5*323);
end

function testVolumeKnownGap(tc)
v = speed_reducer_analysis(xval(), speed_reducer_params('validation2400')).display.V_cm3;
verifyEqual(tc, v, 3958.95, 'AbsTol', 0.05);        % brief's formula
verifyGreaterThan(tc, abs(v - 4147)/4147, 0.04);    % documents the unresolved 4.5 % gap
end

function testBadSize(tc)
verifyError(tc, @() speed_reducer_analysis(ones(6,1), speed_reducer_params()), 'speed_reducer_analysis:badSize');
end

function testNaN(tc)
x = xval();  x(1) = NaN;
verifyError(tc, @() speed_reducer_analysis(x, speed_reducer_params()), 'speed_reducer_analysis:nonFinite');
end

function testNonPositive(tc)
x = xval();  x(3) = 0;
verifyError(tc, @() speed_reducer_analysis(x, speed_reducer_params()), 'speed_reducer_analysis:nonPositive');
end

function testBadCase(tc)
verifyError(tc, @() speed_reducer_params('nonsense'), 'speed_reducer_params:badCase');
end

function testStudentIdParameters(tc)
p = speed_reducer_params('optimisation');
verifyEqual(tc, p.T, 2400);
verifyEqual(tc, p.u, 3.2, 'AbsTol', 1e-12);
end

function testScalingRoundTrip(tc)
p = speed_reducer_params();
xs = [0.1; 0.9; 0.5; 0.3; 0.7; 0.2; 0.6];
verifyEqual(tc, speed_reducer_scale(speed_reducer_scale(xs, p, 'from'), p, 'to'), xs, 'AbsTol', 1e-14);
end

function testComplexStepMatchesCentralDifference(tc)
p = speed_reducer_params();
xs = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7];
[~, ~, gc] = speed_reducer_constraints(xs, p);
h = 1e-6;  Jd = zeros(11, 7);
for k = 1:7
    e = zeros(7,1);  e(k) = h;
    Jd(:,k) = (speed_reducer_constraints(xs+e, p) - speed_reducer_constraints(xs-e, p))/(2*h);
end
verifyEqual(tc, gc.', Jd, 'RelTol', 1e-5, 'AbsTol', 1e-7);
end

function testKnownOptimumFeasible(tc)
p = speed_reducer_params();
x = [3.5; 7.0; 17; 7.3; 7.4; 3.44361; 5.0];
c = speed_reducer_constraints(speed_reducer_scale(x, p, 'to'), p);
verifyLessThan(tc, max(c), 1e-4);
end
