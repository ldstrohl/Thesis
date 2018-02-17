clear all
close all
% load and read data
load PRO2.dat
ca.x = PRO2(:,1);
ca.y = PRO2(:,2);
ca.vx = PRO2(:,3);
ca.vy = PRO2(:,4);
ca.fx = PRO2(:,5);
ca.fy = PRO2(:,6);
ca.fm = PRO2(:,7);
n = length(ca.x);
j=0;

% read background and chaser images
fig = figure;
backg1 = imread('spacebg1.jpg');
backg2 = imread('spacebg2.jpg');
backg3 = imread('spacebg3.jpg');
image(backg1);
chaser1 = imread('chaser1.jpg');
chaser2 = imread('chaser2.jpg');
chaser3 = imread('chaser3.jpg');
fire1 = imread('fire1.jpg');
fire2 = imread('fire2.jpg');
hold on
set(gca,'Visible','off')
fps = 10;
outfile = 'rend.avi';
% aviobj = avifile(outfile,'FPS',fps,'quality',100,'compression','none');
aviobj = VideoWriter(outfile,'Uncompressed AVI');

% rescale the trajectory
ca.xs = ca.x*1.5+870;
ca.ys = ca.y*1.5+510;
for k=1:n
    
    
    % plot chaser and fire
    ca.theta(k) = atan2(-ca.vy(k),ca.vx(k))*180/pi;
    ca.fo(k) = atan2(-ca.fy(k),ca.fx(k))*180/pi;
    if k == n
        ca.theta(k) = 180;
    end
    
    if k<170
        % plot restricted zone and trajectory
        plot(200*cos(0:0.01:2*pi)*1.5/11.6+1115,100*sin(0:0.01:2*pi)*1.5/11.6+475,':','linewidth',1)
        plot(ca.x(1:end)*1.5/11.6+1115,ca.y(1:end)*1.5/11.6+475,'r--');        
        chaser_ro = imrotate(chaser1,ca.theta(k));
        image(ca.x(k)*1.5/11.6+1115,ca.y(k)*1.5/11.6+475-12,chaser_ro);
        %keyboard
%         if k>20
%         keyboard
%         end
%     elseif k>=162 && k<=188
%         % plot restricted zone and trajectory
%         plot(200*cos(0:0.01:2*pi)*1.5+870,100*sin(0:0.01:2*pi)*1.5+510,':','linewidth',1)
%         plot(ca.xs(166:end),ca.ys(166:end),'r--');        
%         chaser_ro = imrotate(chaser2,ca.theta(k));
%         image(ca.xs(k),ca.ys(k)-38,chaser_ro);
    elseif k>170
        % plot restricted zone and trajectory
        plot(400*cos(0:0.01:2*pi)*1.5+870,200*sin(0:0.01:2*pi)*1.5+510,':','linewidth',1)
        plot(ca.xs(171:end)*2-870,ca.ys(171:end)*2-500,'r--');  
        if k>=215
            j=j+1;
            theta_f = atan2(ca.y(215),ca.x(215))*180/pi+180;
            ca.thata(k) = (180-theta_f)/(n-204)*j+theta_f;
            chaser_ro = imrotate(chaser3,ca.thata(k));
        else
            chaser_ro = imrotate(chaser3,ca.theta(k));
        end
        %keyboard
        image(ca.xs(k)*2-870-30,ca.ys(k)*2-38-510,chaser_ro);
        %keyboard
        if ca.fm(k)>=0.2&&ca.fm(k)<=10
            fire1_ro = imrotate(fire1,ca.fo(k));
            if k<=191 %ca.theta(k)<=90 && ca.theta(k)>=0
                image(ca.xs(k)*2-870+(0*cos(ca.theta(k)/180*pi)+50*sin(ca.theta(k)/180*pi)),ca.ys(k)*2-510-0+(-0*sin(ca.theta(k)/180*pi)+50*cos(ca.theta(k)/180*pi)),fire1_ro);
                %keyboard
            elseif k>191 && k<201
                image(ca.xs(k)*2-870-25,ca.ys(k)*2-510+38,fire1_ro);
            else
                image(ca.xs(k)*2-870,ca.ys(k)*2-510+42,fire1_ro);
                %image(ca.xs(k)*2-870+(80*cos(ca.theta(k)/180*pi)-20*sin(ca.theta(k)/180*pi)),ca.ys(k)*2-510+(80*sin(ca.theta(k)/180*pi)+20*cos(ca.theta(k)/180*pi)),fire1_ro);
                %image(ca.xs(k)*2-870+(20*cos(-ca.theta(k)/180*pi)+50*sin(-ca.theta(k)/180*pi)),ca.ys(k)*2-510+(20*sin(-ca.theta(k)/180*pi)+50*cos(-ca.theta(k)/180*pi)),fire1_ro);
               % keyboard
            end
           % keyboard
        elseif ca.fm(k)>=10
            fire2_ro = imrotate(fire2,ca.fo(k));
                image(ca.xs(k)*2-870+5+(0*cos(ca.theta(k)/180*pi)+80*sin(ca.theta(k)/180*pi)),ca.ys(k)*2-510-38+(-0*sin(ca.theta(k)/180*pi)+80*cos(ca.theta(k)/180*pi)),fire2_ro);
                %keyboard
            
        end
    end
 
    hold off
    drawnow;
    frame = getframe(gca);
    aviobj = addframe(aviobj,frame);
    if k<length(ca.xs)
        if k<170
            image(backg1);
%         elseif k>=162 && k<188
%             image(backg2);
        elseif k>=170
            image(backg3);
        end
        set(gca,'Visible','off')
        hold on
    end
end
warning on
aviobj = close(aviobj);