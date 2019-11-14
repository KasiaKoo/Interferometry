classdef PSI
    properties
        al;
        t;
        r;
        wrap;
        notilt;
        phase;
        c_nm;
        Int;
    end
    
    methods
        function obj = PSI(test,reference,algorithm,g,how,n_p,epsilon_p,n_I,epsilon_I)
            if (nargin>0)
                obj.t = test;
                obj.r = reference;
                obj.al = algorithm;
                if algorithm == 3
                    obj=obj.threestep(g,how,n_p,epsilon_p,n_I,epsilon_I);
                elseif algorithm ==4
                    obj = obj.fourstep(g,how,n_p,epsilon_p,n_I,epsilon_I);
                elseif algorithm ==1
                    obj = obj.carre(g,how,n_p,epsilon_p,n_I,epsilon_I);
                elseif algorithm ==2
                    obj = obj.aver(g,how,n_p,epsilon_p,n_I,epsilon_I);                    
                end
                
            end               
           
        end
        
        function obj = threestep(obj,g,how,n_p,epsilon_p,n_I,epsilon_I)
            
            alpha_des = pi/4;
            alpha = phase_error(alpha_des, n_p, epsilon_p);
            I2=abs(exp(1i*3*alpha).*obj.r.front+obj.t.front).^2;
            I2 = detector_error(I2,n_I,epsilon_I);
            
            alpha = phase_error(alpha_des, n_p, epsilon_p);
            I3=abs(exp(1i*5*alpha).*obj.r.front + obj.t.front).^2;
            I3 = detector_error(I3,n_I,epsilon_I);
            
            alpha = phase_error(alpha_des, n_p, epsilon_p);
            I1=abs(exp(1i*(alpha)).*obj.r.front + obj.t.front).^2;
            I1 = detector_error(I1,n_I,epsilon_I);
            
            obj.wrap=g.circa.*(atan2((I3-I2),(I1-I2)));
            %subplot(1,3,3)
            %figure('Name','Wrapped Phase 3','NumberTitle','off')
            %imshow(obj.wrap,[])
            %title('Wrapped Phase 3')
            obj = obj.unwrap_phase(how,g);
            obj.Int = I2;
        end
            
        function obj = fourstep(obj,g,how,n_p,epsilon_p,n_I,epsilon_I)
            
            alpha_des = pi/2;
            delta = phase_error(alpha_des, n_p, epsilon_p);
            I1=abs(exp(1i*(0)*delta).*obj.r.front+obj.t.front).^2;
            I1 = detector_error(I1,n_I,epsilon_I);
            
            delta = phase_error(alpha_des, n_p, epsilon_p);
            I2=abs((exp(1i*(1)*delta).*obj.r.front + obj.t.front)).^2;
            I2 = detector_error(I2,n_I,epsilon_I);
            
            delta = phase_error(alpha_des, n_p, epsilon_p);            
            I3=abs(exp(1i*2*delta).*obj.r.front + obj.t.front).^2;
            I3 = detector_error(I3,n_I,epsilon_I);
            
            delta = phase_error(alpha_des, n_p, epsilon_p);            
            I4=abs(exp(1i*3*delta).*obj.r.front + obj.t.front).^2;
            I4 = detector_error(I4,n_I,epsilon_I);
            
            obj.wrap=(1).*g.circa.*atan2((I4-I2),(I1-I3));
            %subplot(1,3,3)
            %figure('Name','Wrapped Phase 4','NumberTitle','off')
            %imshow(obj.wrap,[])
            %title('Wrapped Phase 4')
            obj.Int = I1;
            obj = obj.unwrap_phase(how,g);
        end
        

        function obj = carre(obj,g,how,n_p,epsilon_p,n_I,epsilon_I)
            
            alpha_des = pi/8;
            
            delta = phase_error(alpha_des, n_p, epsilon_p);
            I1=abs((exp(1i*(-3*delta)).*obj.r.front+obj.t.front)).^2;
            I1 = detector_error(I1,n_I,epsilon_I);
            
            delta = phase_error(alpha_des, n_p, epsilon_p);
            I2=abs((exp(1i*(-delta)).*obj.r.front + obj.t.front)).^2;
            I2 = detector_error(I2,n_I,epsilon_I);
            
            delta = phase_error(alpha_des, n_p, epsilon_p);            
            I3=abs(exp(1i*delta).*obj.r.front + obj.t.front).^2;
            I3 = detector_error(I3,n_I,epsilon_I);
            
            delta = phase_error(alpha_des, n_p, epsilon_p);            
            I4=abs(exp(1i*3*delta).*obj.r.front + obj.t.front).^2;
            I4 = detector_error(I4,n_I,epsilon_I);
            
            top = sqrt(abs((3.*(I2-I3)-(I1-I4)).*((I1-I4)+(I2-I3))));
            bottom = (I2+I3) - (I1+I4);
            obj.wrap=g.circa.*atan2(top,bottom);
            obj.Int = I1;
            obj = obj.unwrap_phase(how,g);
        end
        
        function obj = aver(obj,g,how,n_p,epsilon_p,n_I,epsilon_I)
            
            alpha_des = pi/4;
            delta = phase_error(alpha_des, n_p, epsilon_p);
            I1=abs(exp(1i*(-3*delta).*obj.r.front+obj.t.front)).^2;
            I1 = detector_error(I1,n_I,epsilon_I);
            
            delta = phase_error(alpha_des, n_p, epsilon_p);
            I2=abs((exp(1i*(-1*delta)).*obj.r.front + obj.t.front)).^2;
            I2 = detector_error(I2,n_I,epsilon_I);
            
            delta = phase_error(alpha_des, n_p, epsilon_p);            
            I3=abs(exp(1i*(delta)).*obj.r.front + obj.t.front).^2;
            I3 = detector_error(I3,n_I,epsilon_I);
            
            delta = phase_error(alpha_des, n_p, epsilon_p);            
            I4=abs(exp(1i*3*delta).*obj.r.front + obj.t.front).^2;
            I4 = detector_error(I4,n_I,epsilon_I);
            
            t1 = atan2((1-cos(alpha_des)).*(I1-I3),(sin(alpha_des)).*(2*I3-I1-I3));
            t2 = atan2((1-cos(alpha_des)).*(I2-I4),(sin(alpha_des)).*(2*I3-I2-I4));
            obj.wrap=g.circa.*((t1+t2)./2);
            obj.Int = I1;
            obj = obj.unwrap_phase(how,g);
        end        
        
        
        function obj = unwrap_phase(obj, how,g)
            if how == 'centre'
                n = 1/g.resolution;
                q=obj.wrap;
                q(n:end,:)=unwrap(q(n:end,:),[],1);
                q=flipud(q);
                q(n:end,:)=unwrap(q(n:end,:),[],1);
                q(:,n:end)=unwrap(q(:,n:end),[],2);
                q=fliplr(q);
                q(:,n:end)=unwrap(q(:,n:end),[],2);
                obj.phase = rot90(q,2);
                obj.notilt = obj.phase;
                %subplot(1,3,1)
                %figure('Name','Unwrapped Phase','NumberTitle','off')
                %imshow(obj.phase,[]) 
                %title('Unwrapped Phase')
                obj = obj.zernike_coefficients(g);
                obj = obj.ref_takeaway(g);
                %obj = obj.zernike_coefficients(g);
            elseif how =='edge'
                obj.phase = unwrap(unwrap(p,[],2),[],1);
                obj = obj.ref_takeaway(g);
                obj = obj.zernike_coefficients(g);
            end
        end
        
        function obj = ref_takeaway(obj,g)
            obj.phase = g.circa.*(obj.r.wavefront_error - obj.phase);
            %subplot(1,3,2)
            %figure('Name','Tilt Removed','NumberTitle','off')
            %imshow(obj.phase,[])
            %title('Tilt Removed')
               
        end
        
        function obj = zernike_coefficients(obj,g)
            %Z_00 = 1.*g.circa; %Piston c_00 is calculated as the mean
            Z_11 = 2.*g.r.*cos(g.theta).*g.circa;
            Z_1n1 = 2.*g.r.*sin(g.theta).*g.circa;
            Z_20 = sqrt(3).*(2.*(g.r.^2) - 1).*g.circa;
            Z_22 = sqrt(6).*(g.r.^2).*cos(2.*g.theta).*g.circa;
            Z_2n2 = sqrt(6).*(g.r.^2).*sin(2.*g.theta).*g.circa;
            Z_31 = 2*sqrt(2).*(3.*(g.r.^3) - 2.*g.r).*cos(g.theta).*g.circa;
            Z_3n1 = 2*sqrt(2).*(3.*(g.r.^3) - 2.*g.r).*sin(g.theta).*g.circa;
            Z_33 = 2*sqrt(2).*(g.r.^3).*cos(3.*g.theta).*g.circa;
            Z_3n3 = 2*sqrt(2).*(g.r.^3).*sin(3.*g.theta).*g.circa;
            Z_40 = sqrt(5).*(6.*(g.r.^4) - 6.*(g.r.^2)+1).*g.circa;
            c_lst = linspace(0,0,10);
            abb_options = {Z_11, Z_1n1, Z_20, Z_22, Z_2n2, Z_31, Z_3n1, Z_33, Z_3n3, Z_40};
            abb_name = {'Tip','Tilt','Defocus','Astigmatism Positive'...
                ,'Astigmatism Negative','Coma Positive','Coma Negative','Trefoil Positive','Trefoil Negative','Spherical'};
            for i=1:10
                abe = cell2mat(abb_options(i));
                height = obj.phase.*obj.t.lambda./(4*pi);
                c_lst(i) = nanmean(height(:).*abe(:));
            end
            
            temp_map = containers.Map(abb_name,c_lst);
            obj.c_nm = [keys(temp_map);values(temp_map)]';
        end 
        function obj = plot_interference(obj,g)
            xh = -1:g.resolution:1;
            %figure('Name','Pure Interference','NumberTitle','off')
            image(xh,xh,obj.Int,'CDataMapping','scaled')
            set(gca,'visible','off')
            colormap(gray)
            shading interp 
        end
        function obj = plot_phase(obj,g,which)
            circular_dist = obj.phase;
            if which == 1
                figure();
                surfl(g.x,g.y,circular_dist);
                set(gca,'visible','off');
                colormap(gray);
                shading interp ;
            elseif which ==0
                figure();
                imshow(circular_dist,[]);
            elseif which ==2
                subplot(1,2,1)
                surfl(g.x,g.y,circular_dist);
                set(gca,'visible','off');
                colormap(gray);
                shading interp ;
                subplot(1,2,2)
                imshow(circular_dist,[]);
            end
            
            
        end
        
    end
    
end

    




