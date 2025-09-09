classdef B2KAllParametersLiftOff
    % test = B2KCalc.AllParameters(DimVar(UC(1.131,0.003), 'mm'),DimVar(UC(1.697,0.003), 'mm'),'Water')
    % test = B2KCalc.AllParameters(DimVar(UC(1.131,0.003), 'mm'),DimVar(UC(1.697,0.003), 'mm'),'Water',EXP_flowrate)
    % test = B2KCalc.AllParameters(DimVar(UC(1.131,0.003), 'mm'),DimVar(UC(1.697,0.003), 'mm'),EXP_solvent,EXP_flowrate)
    % test = B2KCalc.AllParameters(DimVar(UC(1.131,B2KCalc.AllParameters.ChannelHeightUncertainty), 'mm'),DimVar(UC(1.697,B2KCalc.AllParameters.ChannelWidthUncertainty), 'mm'),EXP_solvent,EXP_flowrate)
    % test = B2KCalc.AllParameters(DimVar(UC(1.131,B2KCalc.AllParameters.ChannelHeightUncertainty), 'mm'),DimVar(UC(1.697,B2KCalc.AllParameters.ChannelWidthUncertainty), 'mm'),"Water",15)
    properties
        ChannelHeight
        ChannelWidth
        Solvent
        CritFlowRate
    end

    % properties (Constant) 
    % %taken care of in +B2KConstants
    % end

    properties (Constant, Hidden = true) 
        ChannelHeightUncertainty = 0.03; %[mm]
        ChannelWidthUncertainty = 0.03; %[mm]
        CritFlowRateUncertainty = 0.5; %[mL/min]
        % GravitationalAccelerationConstant = DimVar(9.81,'m/s^2')
        GravitationalAccelerationConstant = DimVar(UC(9.81,0),'m/s^2')
        % ChannelLength = DimVar(UC(200,0.003),'mm');
    end

    properties (Constant)
        % ChannelHeightUncertainty = 0.003; %[mm]
        % ChannelWidthUncertainty = 0.003; %[mm]
        % GravitationalAccelerationConstant = DimVar(9.81,'m/s^2')
        ChannelLength = DimVar(UC(200,0.003),'mm');
        RotationalDiameter = DimVar(UC(B2KConstants.pChip.D_max.Value,0),'m'); %[mm]
        EquivVolSpherDiameter = DimVar(UC(B2KConstants.pChip.EquivVolSpherD.Value,0),'m'); %[mm]
    end

    properties (Dependent)
        HydraulicDiameter
        ChannelAspectRatio
        RatioD_maxChannelWidth
        RatioChannelWidthD_max
        RatioD_maxChannelHeight
        RatioChannelHeightD_max
        RatioD_maxSqandChannelHeightChannelWidth
        RatioChannelHeightChannelWidthandD_maxSq
        RatioD_VChannelWidth
        RatioChannelWidthD_V
        RatioD_VChannelHeight
        RatioChannelHeightD_V
        RatioD_VSqandChannelHeightChannelWidth
        RatioChannelHeightChannelWidthandD_VSq
        GravityForce
        BuoyantForce
        ArchimedesNumberUsingD_max
        ArchimedesNumberUsingD_V
        AverageFluidVelocityCRIT
        MaximumFluidVelocityCRIT
        ReynoldsNumberUsingU_avgCRIT
        ParticleReynoldsNumberUsingU_avgD_maxCRIT
        ParticleReynoldsNumberUsingU_avgD_VCRIT
        EntranceLengthTruskeyRectangularDuct
        EntranceLengthBergmanPipe
        TimeToSteadyState
        AvgWallShearRateCRIT
        ShearReynoldsNumberUsingD_maxCRIT
        ShearReynoldsNumberUsingD_VCRIT
        ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT
        ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT
        ModifiedArchimedesNumberUsingD_V
        FroudeNumberUsingU_avgD_maxCRIT
    end
    
    methods
        %--Constructor to initialize object and its base properties--
        function obj = B2KAllParametersLiftOff(channelHeight,channelWidth,solvent,critFlowRate)
            % ChannelHeight
            if matches(channelHeight.UnitStr,'mm')
                channelHeight.Convert('m')
            end
            obj.ChannelHeight = channelHeight;
            
            % ChannelWidth
            if matches(channelWidth.UnitStr,'mm')
                channelWidth.Convert('m')
            end
            obj.ChannelWidth = channelWidth;

            % Solvent
            obj.Solvent = solvent;

            % % FluidFlowRate
            % obj.FluidFlowRateCRIT = DimVar(UC(fluidFlowRate,fluidFlowRate*0.001),'mL/min');
            % %need to add error for fluid flow rate estimate for critical flow rate

            % FluidFlowRate
            obj.CritFlowRate = critFlowRate;
        end

        % %[DO NOT USE SET METHOD TO VALIDATE PROPERTY, ALREADY HANDLED IN DIMVAR, UC]
        % %--Set method to error check set of ChannelHeight
        % function obj = set.ChannelHeight(obj,val)
        %     if val < 0 
        %         error('Channel Height must be positive')
        %     end
        %     obj.ChannelHeight = val;
        % end
        % 
        % %--Set method to error check set of channelheight
        % function obj = set.ChannelWidth(obj,val)
        %     if val < 0 
        %         error('Channel Width must be positive')
        %     end
        %     obj.ChannelWidth = val;
        % end

        %--Get Hydraulic Diameter--
        function val = get.HydraulicDiameter(obj)
            val = 2 .* (obj.ChannelHeight .* obj.ChannelWidth) ./ (obj.ChannelHeight + obj.ChannelWidth);
        end

        %--Get Channel Aspect Ratio--
        function val = get.ChannelAspectRatio(obj)
            val = obj.ChannelHeight ./ obj.ChannelWidth;
        end      

        %--Get RatioD_maxChannelWidth--
        function val = get.RatioD_maxChannelWidth(obj)
            val = B2KConstants.pChip.D_max ./ obj.ChannelWidth;
        end

        %--Get RatioChannelWidthD_max--
        function val = get.RatioChannelWidthD_max(obj)
            val = obj.ChannelWidth ./ B2KConstants.pChip.D_max;
        end

        %--Get RatioD_maxChannelHeight--
        function val = get.RatioD_maxChannelHeight(obj)
            val = B2KConstants.pChip.D_max ./ obj.ChannelHeight;
        end

        %--Get RatioChannelHeightD_max--
        function val = get.RatioChannelHeightD_max(obj)
            val = obj.ChannelHeight ./ B2KConstants.pChip.D_max;
        end

        %--Get RatioD_maxSqandChannelHeightChannelWidth--
        function val = get.RatioD_maxSqandChannelHeightChannelWidth(obj)
            val = B2KConstants.pChip.D_max.^2 ./ (obj.ChannelHeight .* obj.ChannelWidth);
        end

        %--Get RatioChannelHeightChannelWidthandD_maxSq--
        function val = get.RatioChannelHeightChannelWidthandD_maxSq(obj)
            val = obj.ChannelHeight .* obj.ChannelWidth ./ B2KConstants.pChip.D_max.^2;
        end

        %--Get RatioD_VChannelWidth--
        function val = get.RatioD_VChannelWidth(obj)
            val = B2KConstants.pChip.EquivVolSpherD ./ obj.ChannelWidth;
        end

        %--Get RatioChannelWidthD_V--
        function val = get.RatioChannelWidthD_V(obj)
            val = obj.ChannelWidth ./ B2KConstants.pChip.EquivVolSpherD;
        end

        %--Get RatioD_VChannelHeight--
        function val = get.RatioD_VChannelHeight(obj)
            val = B2KConstants.pChip.EquivVolSpherD ./ obj.ChannelHeight;
        end

        %--Get RatioChannelHeightD_V--
        function val = get.RatioChannelHeightD_V(obj)
            val = obj.ChannelHeight ./ B2KConstants.pChip.EquivVolSpherD;
        end

        %--Get RatioD_VSqandChannelHeightChannelWidth--
        function val = get.RatioD_VSqandChannelHeightChannelWidth(obj)
            val = B2KConstants.pChip.EquivVolSpherD.^2 ./ (obj.ChannelHeight .* obj.ChannelWidth);
        end

        %--Get RatioChannelHeightChannelWidthandD_VSq--
        function val = get.RatioChannelHeightChannelWidthandD_VSq(obj)
            val = obj.ChannelHeight .* obj.ChannelWidth ./ B2KConstants.pChip.EquivVolSpherD.^2;
        end

        %--Get GravitationalForce--
        function val = get.GravityForce(obj)
            % val = B2KConstants.pChip.Density .* (B2KConstants.pChip.H .* B2KConstants.pChip.W .* B2KConstants.pChip.W) .* (obj.GravitationalAccelerationConstant);
            val = B2KConstants.pChip.Density .* (B2KConstants.pChip.H .* B2KConstants.pChip.W .* B2KConstants.pChip.W) .* (obj.GravitationalAccelerationConstant);
        end

        %--Get BuoyantForce--
        function val = get.BuoyantForce(obj)
            val = B2KConstants.(obj.Solvent).rho .* (B2KConstants.pChip.H .* B2KConstants.pChip.W .* B2KConstants.pChip.W) .* (obj.GravitationalAccelerationConstant);
        end

        %--Get ArchimedesNumberUsingD_max--
        function val = get.ArchimedesNumberUsingD_max(obj)
            val = (B2KConstants.(obj.Solvent).rho .* (B2KConstants.pChip.Density - B2KConstants.(obj.Solvent).rho) .* (obj.GravitationalAccelerationConstant) .* (B2KConstants.pChip.D_max).^3) ./ (B2KConstants.(obj.Solvent).mu).^2;
        end

        %--Get ArchimedesNumberUsingEquivVolSpherD--
        function val = get.ArchimedesNumberUsingD_V(obj)
            val = (B2KConstants.(obj.Solvent).rho .* (B2KConstants.pChip.Density - B2KConstants.(obj.Solvent).rho) .* (obj.GravitationalAccelerationConstant) .* (B2KConstants.pChip.EquivVolSpherD).^3) ./ (B2KConstants.(obj.Solvent).mu).^2;
        end

        %--Get AverageFluidVelocityCRIT--
        function val = get.AverageFluidVelocityCRIT(obj)
            val = obj.CritFlowRate ./ (obj.ChannelHeight .* obj.ChannelWidth);
        end

        %--Get MaximumFluidVelocityCRIT--
        function val = get.MaximumFluidVelocityCRIT(obj)
            tol = 1e-10;

            % Initialize sums and error estimates
            sum = 0;
            error = 1;
            n = 0;

            % Loop to calculate the infinite series with dynamic termination
            while error > tol
                term = tanh((2 .* n + 1) .* pi .* obj.ChannelWidth.Value.Value ./ (2 .* obj.ChannelHeight.Value.Value)) ./ ((2 .* n + 1).^5 * pi.^5);
                sum = sum + term;
                error = abs(term); % Update error estimate
                n = n + 1;
            end
        
            % Calculate the pressure drop
            denominator = obj.ChannelHeight.^3 .* obj.ChannelWidth .* (1 - 6 .* (obj.ChannelHeight ./ obj.ChannelWidth) .* sum);
            pressuredrop_div_L = (12 .* B2KConstants.(obj.Solvent).mu .* obj.CritFlowRate) ./ denominator;

            val = (obj.ChannelHeight.^2 ./ (8 .* B2KConstants.(obj.Solvent).mu)) .* pressuredrop_div_L;
        end

        %--Get ReynoldsNumberUsingU_avgCRIT--
        function val = get.ReynoldsNumberUsingU_avgCRIT(obj)
            val = B2KConstants.(obj.Solvent).rho .* obj.AverageFluidVelocityCRIT .* obj.HydraulicDiameter ./ B2KConstants.(obj.Solvent).mu;
        end

        %--Get ParticleReynoldsNumberUsingU_avgD_maxCRIT--
        function val = get.ParticleReynoldsNumberUsingU_avgD_maxCRIT(obj)
            val = obj.ReynoldsNumberUsingU_avgCRIT .* (B2KConstants.pChip.D_max ./ obj.HydraulicDiameter).^2;
        end

        %--Get ParticleReynoldsNumberUsingU_avgD_VCRIT--
        function val = get.ParticleReynoldsNumberUsingU_avgD_VCRIT(obj)
            val = obj.ReynoldsNumberUsingU_avgCRIT .* (B2KConstants.pChip.EquivVolSpherD ./ obj.HydraulicDiameter).^2;
        end

        %--Get EntranceLengthTruskeyRectangularDuct--
        function val = get.EntranceLengthTruskeyRectangularDuct(obj)
            val = 0.04 .* obj.ChannelHeight .* obj.ReynoldsNumberUsingU_avgCRIT;
        end

        %--Get EntranceLengthBergmanPipe--
        function val = get.EntranceLengthBergmanPipe(obj)
            val = 0.058 .* obj.HydraulicDiameter .* obj.ReynoldsNumberUsingU_avgCRIT;
        end

        %--Get TimeToSteadyState--
        function val = get.TimeToSteadyState(obj)
            val = obj.EntranceLengthTruskeyRectangularDuct ./ obj.AverageFluidVelocityCRIT;
        end

        %--Get AvgWallShearRateCRIT--
        function val = get.AvgWallShearRateCRIT(obj)
            % % syms n
            % % gamma1 = (12 .* B2KConstants.(obj.Solvent).mu .* obj.ChannelLength .* obj.FluidFlowRate) ./ (obj.ChannelHeight .^3 .* obj.ChannelWidth .*(1 - 6 .* (obj.ChannelHeight ./ obj.ChannelWidth) .* symsum((tanh((2 .* n + 1) .* pi .* obj.ChannelWidth ./ (2 .* obj.ChannelHeight))) ./ ((2 .* n + 1) .^5 .* pi .^5),n,[0 Inf])));
            % % gamma2 = (obj.ChannelHeight ./ (2 .* B2KConstants.(obj.Solvent).mu .* obj.ChannelLength)) .* (1 - 16 .* (obj.ChannelHeight ./ obj.ChannelWidth) .* symsum((tanh((2 .* n + 1) .* pi .* obj.ChannelWidth ./ (2 .* obj.ChannelHeight))) ./ ((2 .* n + 1) .^5 .* pi .^5),n,[0 Inf]));
            % % val = double(gamma1 .* gamma2);
            % % val = double((12 .* B2KConstants.(obj.Solvent).mu .* obj.ChannelLength .* obj.FluidFlowRate) ./ (obj.ChannelHeight .^3 .* obj.ChannelWidth .*(1 - 6 .* (obj.ChannelHeight ./ obj.ChannelWidth) .* symsum((tanh((2 .* n + 1) .* pi .* obj.ChannelWidth ./ (2 .* obj.ChannelHeight))) ./ ((2 .* n + 1) .^5 .* pi .^5),n,[0 Inf]))) .* (obj.ChannelHeight ./ (2 .* B2KConstants.(obj.Solvent).mu .* obj.ChannelLength)) .* (1 - 16 .* (obj.ChannelHeight ./ obj.ChannelWidth) .* symsum((tanh((2 .* n + 1) .* pi .* obj.ChannelWidth ./ (2 .* obj.ChannelHeight))) ./ ((2 .* n + 1) .^5 .* pi .^5),n,[0 Inf])));
            % 
            % % n = sym('n');
            % % gamma1 = (12 .* B2KConstants.(obj.Solvent).mu .* obj.ChannelLength .* obj.FluidFlowRate) ./ (obj.ChannelHeight .^3 .* obj.ChannelWidth .*(1 - 6 .* (obj.ChannelHeight ./ obj.ChannelWidth) .* symsum((tanh((2 .* n + 1) .* pi .* obj.ChannelWidth ./ (2 .* obj.ChannelHeight))) ./ ((2 .* n + 1) .^5 .* pi .^5),n,[0 Inf])));
            % % gamma2 = (obj.ChannelHeight ./ (2 .* B2KConstants.(obj.Solvent).mu .* obj.ChannelLength)) .* (1 - 16 .* (obj.ChannelHeight ./ obj.ChannelWidth) .* symsum((tanh((2 .* n + 1) .* pi .* obj.ChannelWidth ./ (2 .* obj.ChannelHeight))) ./ ((2 .* n + 1) .^5 .* pi .^5),n,[0 Inf]));
            % % val = double(gamma1 .* gamma2);
            % 
            % %**symbolic toolbox probably can't accomodate DimVar and UC**
            % 
            % %--Using approximation for AvgWallShearRateCRIT without use of Symbolic Toolbox
            % val = 6 .* obj.FluidFlowRate ./ (obj.ChannelWidth .* obj.ChannelHeight .^2);

            % %
            tol = 1e-10;

            % Initialize sums and error estimates
            sum1 = 0;
            sum2 = 0;
            error1 = 1;
            error2 = 1;
            n = 0;

            % Loop to calculate the infinite series with dynamic termination
            while error1 > tol || error2 > tol
                % term1 = tanh((2 .* n + 1) .* (pi .* obj.ChannelWidth) / (2 .* obj.ChannelHeight)) ./ (((2 .* n + 1) .^ 5) .* (pi .^ 5));
                term1 = tanh((2 .* n + 1) .* (pi .* obj.ChannelWidth.Value.Value) / (2 .* obj.ChannelHeight.Value.Value)) ./ (((2 .* n + 1) .^ 5) .* (pi .^ 5));
                sum1 = sum1 + term1;
                error1 = abs(term1); % Update error estimate for sum1

                % term2 = ((-1) .* tanh((2 .* n + 1) .* (pi .* obj.ChannelWidth) / (2 .* obj.ChannelHeight))) ./ ((2 .* n + 1) .^ 3) .* (pi .^ 3);
                term2 = ((-1) .* tanh((2 .* n + 1) .* (pi .* obj.ChannelWidth.Value.Value) / (2 .* obj.ChannelHeight.Value.Value))) ./ ((2 .* n + 1) .^ 3 .* pi .^ 3);
                sum2 = sum2 + term2;
                error2 = abs(term2); % Update error estimate for sum2

                n = n + 1;
            end

            % Calculate the first part of the equation
            part1 = (6 .* obj.CritFlowRate) ./ (obj.ChannelHeight .^ 2 .* obj.ChannelWidth .* (1 - 6 .* (obj.ChannelHeight ./ obj.ChannelWidth) .* sum1));

            % Calculate the second part of the equation
            part2 = (1 - 16 .* (obj.ChannelHeight ./ obj.ChannelWidth) .* sum2);

            % Calculate gamma_x
            val = part1 .* part2;
        end

        %--Get ShearReynoldsNumberUsingD_max--
        function val = get.ShearReynoldsNumberUsingD_maxCRIT(obj)
            val = (B2KConstants.(obj.Solvent).rho .* obj.AvgWallShearRateCRIT .* B2KConstants.pChip.D_max .^2) ./ B2KConstants.(obj.Solvent).mu;
        end

         %--Get ShearReynoldsNumberUsingD_max--
        function val = get.ShearReynoldsNumberUsingD_VCRIT(obj)
            val = (B2KConstants.(obj.Solvent).rho .* obj.AvgWallShearRateCRIT .* B2KConstants.pChip.EquivVolSpherD .^2) ./ B2KConstants.(obj.Solvent).mu;
        end

        %--Get ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT--
        function val = get.ModifiedParticleReynoldsNumberUsingU_avgD_maxCRIT(obj)
            % num = B2KConstants.(obj.Solvent).rho .* obj.AverageFluidVelocityCRIT .* obj.RotationalDiameter;
            % denom = B2KConstants.(obj.Solvent).mu .* (1.4 - 0.8 .* exp((obj.HydraulicDiameter ./ DimVar(UC(50,0.003),'mm'))/1.5));
            % val = num ./ denom;
            num = obj.ParticleReynoldsNumberUsingU_avgD_maxCRIT;
            denom = (1.25 - 0.5 .* exp((obj.HydraulicDiameter ./ DimVar(UC(50,0.003),'mm'))/1.5));
            val = num ./ denom;
        end

        %--Get ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT--
        function val = get.ModifiedParticleReynoldsNumberUsingU_avgD_VCRIT(obj)
            % num = B2KConstants.(obj.Solvent).rho .* obj.AverageFluidVelocityCRIT .* B2KConstants.pChip.EquivVolSpherD;
            % denom = B2KConstants.(obj.Solvent).mu .* (1.4 - 0.8 .* exp((obj.HydraulicDiameter ./ DimVar(UC(50,0.003),'mm'))/1.5));
            % val = num ./ denom;
            num = obj.ParticleReynoldsNumberUsingU_avgD_VCRIT;
            denom = (1.25 - 0.5 .* exp(-(obj.HydraulicDiameter ./ DimVar(UC(50,0.003),'mm'))/1.5));
            val = num ./ denom;
        end

        %--Get ModifiedArchimedesNumberUsingD_V--
        function val = get.ModifiedArchimedesNumberUsingD_V(obj)
            val = obj.ArchimedesNumberUsingD_V .* 0.70;
        end

        %--Get FroudeNumberUsingU_avgD_maxCRIT--
        function val = get.FroudeNumberUsingU_avgD_maxCRIT(obj)
            val = obj.AverageFluidVelocityCRIT ./ sqrt(obj.GravitationalAccelerationConstant .* B2KConstants.pChip.D_max);
        end
    end
end