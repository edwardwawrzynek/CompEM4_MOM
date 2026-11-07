classdef MOMProblem
    properties
        k0 % wavenumber of free space
        nodes % 2xn matrix of node locations on cylinder
        Z % system matrix (nxn)
    end

    methods
        % initialize 
        function obj = MOMProblem(k0)
            obj.k0 = k0;
        end
        
        % create uniformly distributed nodes over a circle of radius r
        function obj = makeCircle(obj, r, npts)
            phi = linspace(0,2*pi,npts+1);
            phi = phi(1:end-1);

            obj.nodes = [r*cos(phi); r*sin(phi)];
            obj.Z = zeros(npts);
        end

        % return the midpoint, length, and normal vector of the n-th pulse basis function
        function [midpt, Ln, n] = pulse_basis_func(obj, n)
            Nbasis = size(obj.nodes, 2);
            
            % get locations of endpoint nodes
            x0 = obj.nodes(:,n);
            x1 = obj.nodes(:,mod(n+1, Nbasis));

            midpt = 0.5 * (x0 + x1);
            Ln = norm(x0 - x1);
            % outward unit normal for the circular boundary segment
            tangent = x1 - x0;
            n = [tangent(2); -tangent(1)] / Ln;
        end

        % create the system matrix for EFIE TMz scattering from a PEC cylinder
        % using pulse basis functions and collocation
        function obj = tmz_pec_sys_matrix(obj)
            % number of basis functions
            Nbasis = size(obj.nodes, 2);
            % impedance of space
            eta = 376.730;
            % Euler–Mascheroni constant
            gamma = -psi(1);
            e = exp(1);

            for n = 1:numel(Nbasis)
                [rho_n, Ln, nn] = obj.pulse_basis_func(n);
                for m = 1:numel(Nbasis)
                    [rho_m, Lm, nm] = obj.pulse_basis_func(m);

                    if n == m
                        % handle self interaction through Hankel function at
                        % small argument
                        obj.Z(m,n) = obj.k0 * eta * Ln / 4 * (1 - 2j/pi * log(obj.k0 * gamma * Ln / 4 / e));
                    else
                        % approximate integral with midpoint rule elsewhere
                        obj.Z(m,n) = obj.k0 * eta * Ln / 4 * besselh(0, 2, obj.k0*norm(rho_m - rho_n));
                    end
                end
            end
        end

        % create the system matrix for MFIE TEz scattering from a PEC cylinder
        % using pulse basis functions and collocation
        function obj = tez_pec_sys_matrix(obj)
            % number of basis functions
            Nbasis = size(obj.nodes, 2);

            for n = 1:numel(Nbasis)
                [rho_n, Ln, nn] = obj.pulse_basis_func(n);
                for m = 1:numel(Nbasis)
                    [rho_m, Lm, nm] = obj.pulse_basis_func(m);

                    if n == m
                        % self interaction from MFIE
                        obj.Z(m,n) = -1/2;
                    else
                        % approximate integral with midpoint rule elsewhere
                        obj.Z(m,n) = obj.k0 * Ln / 4j * besselh(1,2,obj.k0*norm(rho_m-rho_n)) * dot(nn, rho_m-rho_n) / norm(rho_m-rho_n); 
                    end
                end
            end
        end
    end
end