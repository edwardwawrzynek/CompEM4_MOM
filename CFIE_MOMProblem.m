classdef CFIE_MOMProblem
    properties
        k0 % wavenumber of free space
        ki % wavenumber of dielectric
        nodes % 2xn matrix of node locations on cylinder
        Z % system matrix (2nx2n)
        V % incidence vector (2nx1)
        I % solution vector (2nx1)
        F0 % excitation magnitude 
        phi_inc % angle of excitation [rad]
    end

    methods
        % initialize 
        function obj = CFIE_MOMProblem(k0)
            obj.k0 = k0;
        end
        
        % create uniformly distributed nodes over a circle of radius r
        function obj = makeCircle(obj, r, npts)
            phi = linspace(0,2*pi,npts+1);
            phi = phi(1:end-1);

            obj.nodes = [r*cos(phi); r*sin(phi)];
            obj.Z = zeros(2*npts);
        end

        % return the midpoint, length, and normal vector of the n-th pulse basis function
        function [midpt, Ln, n] = pulse_basis_func(obj, n)
            Nbasis = size(obj.nodes, 2);
            
            % get locations of endpoint nodes
            x0 = obj.nodes(:,n);
            x1 = obj.nodes(:,mod(n, Nbasis)+1);

            midpt = 0.5 * (x0 + x1);
            Ln = norm(x0 - x1);
            % outward unit normal for the circular boundary segment
            tangent = x1 - x0;
            n = [tangent(2); -tangent(1)] / Ln;
        end

        function plotMesh(obj)
        % PLOTMESH Plots the mesh segments for a MOMProblem object.
        % obj.plotMesh()
        %
        % INPUTS
        % obj       the MOMProblem object for which the mesh segments are 
        %           plotted
        %
        % OUTPUTS
        % 
            
            % close the distance between the first and last node
            segments = [obj.nodes, obj.nodes(:,1)];
            
            hold on;
            plot(segments(1,:), segments(2,:), 'b-', 'LineWidth', 1.5); % plot segments
            plot(obj.nodes(1,:), obj.nodes(2,:), 'ro', 'MarkerFaceColor', 'r'); % plot nodes

        end

        % create the system matrix for the CFIE TMz scattering from a
        % dielectric cylinder using pulse basis funcions and collocation
        % the first half of the unknown vector obj.I is Ez and the second
        % half is Ht
        function obj = tmz_eps_sys_matrix(obj, eps_r)
            obj.ki = sqrt(eps_r) * obj.k0;
            mu_r = 1;
            % number of basis functions
            Nbasis = size(obj.nodes, 2);
            % impedance of space
            eta = 376.730;
            % Euler–Mascheroni constant
            gamma = exp(-psi(1));
            e = exp(1);
            % wavenumbers in free space and dielectric
            k0 = obj.k0;
            ki = k0 * sqrt(eps_r);

            % Z is divided as
            % top left: Ez to Ezinc through dG0/dn
            % top right: Ht to Ezinc through G0
            % bottom left: Ez to 0 through dGi/dn
            % bottom right: Ht to 0 through Gi

            for n = 1:Nbasis
                [rho_n, Ln, nn] = obj.pulse_basis_func(n);
                for m = 1:Nbasis
                    [rho_m, Lm, nm] = obj.pulse_basis_func(m);

                    if n == m
                        % EXTERNAL interactions (10.2.56)
                        obj.Z(m,n) = 1/2;
                        obj.Z(m,n+Nbasis) = k0 * eta * Ln / 4 * (1 - 2j/pi * log(k0 * gamma * Ln / 4 / e));
                        % INTERNAL interactions (10.2.57)
                        obj.Z(m+Nbasis,n) = 1/2;
                        obj.Z(m+Nbasis,n+Nbasis) = -mu_r * ki * eta * Ln / 4 * (1 - 2j/pi * log(ki * gamma * Ln / 4 / e));
                    else
                        % EXTERNAL interactions (10.2.56)
                        obj.Z(m,n) = -k0 * Ln / 4j * besselh(1,2,k0*norm(rho_m-rho_n)) * dot(nn, rho_m-rho_n) / norm(rho_m-rho_n); 
                        obj.Z(m,n+Nbasis) = k0 * eta * Ln / 4 * besselh(0, 2, k0*norm(rho_m - rho_n));
                        % INTERNAL interactions (10.2.57)
                        obj.Z(m+Nbasis,n) = ki * Ln / 4j * besselh(1,2,ki*norm(rho_m-rho_n)) * dot(nn, rho_m-rho_n) / norm(rho_m-rho_n);
                        obj.Z(m+Nbasis,n+Nbasis) = -mu_r * ki * eta * Ln / 4 * besselh(0, 2, ki*norm(rho_m - rho_n));
                    end
                end
            end
        end

        % create the system matrix for the CFIE TEz scattering from a
        % dielectric cylinder using pulse basis funcions and collocation
        % the first half of the unknown vector obj.I is Hz and the second
        % half is Et
        function obj = tez_eps_sys_matrix(obj, eps_r)
            obj.ki = sqrt(eps_r) * obj.k0;
            mu_r = 1;
            % number of basis functions
            Nbasis = size(obj.nodes, 2);
            % admittance of free space
            Y0 = 1/376.730;
            % Euler–Mascheroni constant
            gamma = exp(-psi(1));
            e = exp(1);
            % wavenumbers in free space and dielectric
            k0 = obj.k0;
            ki = k0 * sqrt(eps_r);

            % Z is divided as
            % top left: Hz to Hzinc through dG0/dn
            % top right: Et to Hzinc through G0
            % bottom left: Hz to 0 through dGi/dn
            % bottom right: Et to 0 through Gi

            for n = 1:Nbasis
                [rho_n, Ln, nn] = obj.pulse_basis_func(n);
                for m = 1:Nbasis
                    [rho_m, Lm, nm] = obj.pulse_basis_func(m);

                    if n == m
                        % EXTERNAL interactions (10.2.56)
                        obj.Z(m,n) = 1/2;
                        obj.Z(m,n+Nbasis) = -k0 * Y0 * Ln / 4 * (1 - 2j/pi * log(k0 * gamma * Ln / 4 / e));
                        % INTERNAL interactions (10.2.57)
                        obj.Z(m+Nbasis,n) = 1/2;
                        obj.Z(m+Nbasis,n+Nbasis) = eps_r * ki * Y0 * Ln / 4 * (1 - 2j/pi * log(ki * gamma * Ln / 4 / e));
                    else
                        % EXTERNAL interactions (10.2.56)
                        obj.Z(m,n) = -k0 * Ln / 4j * besselh(1,2,k0*norm(rho_m-rho_n)) * dot(nn, rho_m-rho_n) / norm(rho_m-rho_n); 
                        obj.Z(m,n+Nbasis) = -k0 * Y0 * Ln / 4 * besselh(0, 2, k0*norm(rho_m - rho_n));
                        % INTERNAL interactions (10.2.57)
                        obj.Z(m+Nbasis,n) = ki * Ln / 4j * besselh(1,2,ki*norm(rho_m-rho_n)) * dot(nn, rho_m-rho_n) / norm(rho_m-rho_n);
                        obj.Z(m+Nbasis,n+Nbasis) = eps_r * ki * Y0 * Ln / 4 * besselh(0, 2, ki*norm(rho_m - rho_n));
                    end
                end
            end
        end

        function obj = constructIncidence(obj, phi_inc, F0)
        % CONSTRUCTINCIDENCE constructs the RHS of the MoM matrix using
        % information about the incident field.
        % obj.constructIncidence(phi_inc, F0)
        %
        % INPUTS
        % phi_inc   the angle of incidence as measured from +x in radians
        % F0        the magnitude of the incident field
        %
        % OUTPUTS
        % obj       the object with updated properties related to
        %           excitation and V, the RHS of the MoM matrix problem

            % save excitation parameters
            obj.phi_inc = phi_inc; 
            obj.F0 = F0;

            % k vector
            k_inc = obj.k0 * [cos(phi_inc); sin(phi_inc)];

            obj.V = zeros(2*size(obj.nodes, 2), 1); % initialize
            for m = 1:size(obj.nodes,2) % for each node
                rho_m = obj.pulse_basis_func(m); % location vector
                obj.V(m) = F0 * exp(-1i * dot(k_inc, rho_m)); % incident E-field
            end

        end

        function obj = solveI(obj)
        % SOLVEI Solves for the currents on the mesh elements pertaining to
        % the MoM object obj.
        % obj.solveI( )
        %
        % INPUTS
        %
        % OUTPUTS
        % obj       the object with an updated solution vector for the
        %           currents I.
            obj.I = obj.Z \ obj.V;
        end

        function [Finc, Fscat, Ftot] = evaluateEPSFields(obj, X, Y, pol)
        % EVALUATEEPSFIELDS Solves for fields as a result of the surface
        % currents in I for both TMz and TEz cases.
        % [Finc, Fscat, Ftot] = obj.evaluateEpsFields(X,Y,pol)
        %
        % INPUTS
        % X         the x-coordinate of points for which the field is
        %           evaluated
        % Y         the y-coordinate of points for which the field is
        %           evaluated
        % pol       a value of 1 if TMz polarization, a value of 2 if TEz
        %           polarization
        %
        % OUTPUTS
        % Finc      the analytical incident field at points (X,Y)
        % Fscat     the scattered field at points (X,Y)
        % Ftot      the total field at points (X,Y)

            % pol = 1 corresponds to TMz
            % pol = 2 corresponds to TEz
            if pol ~= 1 && pol ~= 2
                error("Invalid case: set pol = 1 for TMz polarization or pol = 2 for TEz polarization");
            end

            % initialize
            Fscat = zeros(numel(Y), numel(X));
            Finc = zeros(numel(Y), numel(X));
            Ftot = zeros(numel(Y), numel(X));
            
            Nbasis = size(obj.nodes, 2);

            % free space impedance
            Y0 = 1/376.730;
            eps_r = (obj.ki / obj.k0)^2;

            % for each field point
            for ix = 1:numel(X)
                for iy = 1:numel(Y)
                    fieldPoint = [X(ix); Y(iy)]; % calculate target location
                    interiorPoint = norm(fieldPoint) < norm(obj.nodes(:,1)); % check if point is inside dielectric

                    for n = 1:Nbasis

                        [sourcePoint, Ln, nn] = obj.pulse_basis_func(n); % source info
                        r = fieldPoint - sourcePoint; % distance
                        distance = norm(r); % magnitude of distance

                        if distance > 0
                            if pol == 1
                                % construct the contribution from each node
                                if ~interiorPoint
                                    Fscat(iy, ix) = Fscat(iy, ix) - ...
                                        obj.k0 * eta * Ln / 4 * obj.I(n+Nbasis) * besselh(0, 2, obj.k0 * distance) + ...
                                        obj.k0 * Ln / 4j * obj.I(n) * besselh(1, 2, obj.k0 * distance) * dot(nn, r) / distance;
                                else
                                    Fscat(iy, ix) = Fscat(iy, ix) + ...
                                        obj.ki * eta * Ln / 4 * obj.I(n+Nbasis) * besselh(0, 2, obj.ki * distance) - ...
                                        obj.ki * Ln / 4j * obj.I(n) * besselh(1, 2, obj.ki * distance) * dot(nn, r) / distance;
                                end
                            else
                                % construct the contribution from each node
                                % construct the contribution from each node
                                if ~interiorPoint
                                    Fscat(iy, ix) = Fscat(iy, ix) + ...
                                        obj.k0 * Y0 * Ln / 4 * obj.I(n+Nbasis) * besselh(0, 2, obj.k0 * distance) + ...
                                        obj.k0 * Ln / 4j * obj.I(n) * besselh(1, 2, obj.k0 * distance) * dot(nn, r) / distance;
                                else
                                    Fscat(iy, ix) = Fscat(iy, ix) - ...
                                        obj.ki * eps_r * Y0 * Ln / 4 * obj.I(n+Nbasis) * besselh(0, 2, obj.ki * distance) - ...
                                        obj.ki * Ln / 4j * obj.I(n) * besselh(1, 2, obj.ki * distance) * dot(nn, r) / distance;
                                end
                            end
                        end
                    end

                    % analytical incident field
                    Finc(iy, ix) = obj.F0 * exp(-1j * obj.k0 * ...
                        (X(ix)*cos(obj.phi_inc) + Y(iy)*sin(obj.phi_inc)));

                    % total field
                    if ~interiorPoint
                        Ftot(iy, ix) = Finc(iy, ix) + Fscat(iy, ix);
                    else
                        Ftot(iy, ix) = Fscat(iy, ix);
                    end
                end
            end
        end

    end
end