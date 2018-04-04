%%
function edit_components_GUI(DFF,A,C,b,f,Cn,S,options)

% reject = zeros(size(A,2),1);
rejectFile = matfile('tmp.mat','Writable',true);
rejectFile.reject = zeros(size(A,2),1);

% memmaped = isobject(Y);
defoptions = CNMFSetParms;
if nargin < 7 || isempty(options); options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end          % # of columns
% if ~isfield(options,'normalize') || isempty(options.normalize); options.normalize = ones(size(A,1),1); end
%     sn = options.normalize;
if ~isfield(options,'plot_df') || isempty(options.plot_df); options.df = defoptions.plot_df; end
plot_df = options.plot_df;
if ~isfield(options,'make_gif') || isempty(options.make_gif); options.make_gif = defoptions.make_gif; end
make_gif = options.make_gif;
if ~isfield(options,'save_avi') || isempty(options.save_avi); options.save_avi = defoptions.save_avi; end
save_avi = options.save_avi;
if ~isfield(options,'sx') || isempty(options.sx); options.sx = defoptions.sx; end
sx = min([options.sx,floor(d1/2),floor(d2/2)]);
if isfield(options,'name') && ~isempty(options.name)
    name = [options.name,'_components'];
else
    name = [defoptions.name,'_components'];
end
if ~isfield(options,'full_A') || isempty(options.full_A); full_A = defoptions.full_A; else full_A = options.full_A; end
% if ~memmaped; Y = double(Y); end
T = size(C,2);
% mov =Y;
% if ndims(Y) == 3
%     Y = reshape(Y,d1*d2,T);
% end
if nargin < 6 || isempty(Cn)
    Cn = reshape(mean(Y,2),d1,d2);
end
b = double(b);
C = double(C);
f = double(f);
nA = full(sqrt(sum(A.^2))');
[K,~] = size(C);
A = A/spdiags(nA,0,K,K);    % normalize spatial components to unit energy
C = bsxfun(@times,C,nA(:)); %spdiags(nA,0,K,K)*C;

nr = size(A,2);     % number of ROIs
nb = size(f,1);     % number of background components
%nA = full(sum(A.^2))';  % energy of each row
%Y_r = spdiags(nA,0,nr,nr)\(A'*Y- (A'*A)*C - (A'*full(b))*f) + C; 

% AY = mm_fun(A,Y);
% Y_r = (AY- (A'*A)*C - full(A'*double(b))*f) + C;




% if plot_df
%     [~,Df] = extract_DF_F(Y,A,C,[],options,AY);
% else
%     Df = ones(size(A,2)+1,1);
% end

if save_avi
    vidObj = VideoWriter([name,'.avi']);
    set(vidObj,'FrameRate',1);
    open(vidObj);
end
thr = 0.95;
fig = figure('Visible','off');
set(gcf,'KeyPressFcn',@KeyPressFcn);
set(gcf,'Position',2*[300,300,960,480]);
set(gcf,'PaperPosition',2*[300,300,960,480]);
int_x = zeros(nr,2*sx);
int_y = zeros(nr,2*sx);
cm = com(A,d1,d2);



% Create a figure and axes

% ax = axes('Units','DF/F');



% Create slider
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',nr+nb,'Value',1,'SliderStep',[1/(nr+nb-1) 1],...
    'Position', [150 20 800 20],...
    'Callback', @surfzlim);


% Add a text uicontrol to label the slider.
txt = uicontrol('Style','text',...
    'Position',[400 45 120 20],...
    'String','Component');

% sld_t = uicontrol('Style','slider',...
%     'Min',1,'Max',T,'Value',1,'SliderStep',[1/(T-1) 1],...
%     'Position',[1000 20 800 20],...
%     'Callback',@rawDatUpdate);
% txt_t = uicontrol('Style','text',...
%     'Position',[1200 45 120 20],...
%     'String','Frame #');

rej_but = uicontrol('Style','pushbutton',...
    'String','Reject',...
    'Position',[1600 250 200 50], ...
    'Callback',@rejectComponent);

unrej_but = uicontrol('Style','pushbutton',...
    'String','Un-Reject',...
    'Position',[1600 100 200 50],...
    'Callback',@unrejectComponent);

% % % % add controls for arrow keys and hot key for reject and undo reject
% % %  rearrange plots 


% Make figure visble after adding all components
fig.Visible = 'on';


global x_cont
global y_cont
plot_component(1)
% plotRawDat(1)


% This code uses dot notation to set properties.
% Dot notation runs in R2014b and later.
% For R2014a and earlier: set(f,'Visible','on');


    function surfzlim(source,callbackdata)
        i = source.Value;
        plot_component(round(i));
        
    end

    function plot_component(i)
       if i <= nr
            subplot(3,4,9);
            Atemp = reshape(A(:,i),d1,d2);
            int_x(i,:) = round(cm(i,1)) + (-(sx-1):sx);
            if int_x(i,1)<1
                int_x(i,:) = int_x(i,:) + 1 - int_x(i,1);
            end
            if int_x(i,end)>d1
                int_x(i,:) = int_x(i,:) - (int_x(i,end)-d1);
            end
            int_y(i,:) = round(cm(i,2)) + (-(sx-1):sx);
            if int_y(i,1)<1
                int_y(i,:) = int_y(i,:) + 1 - int_y(i,1);
            end
            if int_y(i,end)>d2
                int_y(i,:) = int_y(i,:) - (int_y(i,end)-d2);
            end
            Atemp = Atemp(int_x(i,:),int_y(i,:));
            imagesc(int_x(i,:),int_y(i,:),Atemp); axis square;
        end
        subplot(3,4,[1, 2,5, 6]);
        if i <= nr
            cla
            imagesc(2*Cn); hold on;
            A_temp = full(reshape(A(:,i),d1,d2));
            A_temp = medfilt2(A_temp,[3,3]);
            A_temp = A_temp(:);
            [temp,ind] = sort(A_temp(:).^2,'ascend');
            temp =  cumsum(temp);
            ff = find(temp > (1-thr)*temp(end),1,'first');
            if ~isempty(ff)
                x_cont =  reshape(A_temp,d1,d2);
                y_cont = [0,0]+A_temp(ind(ff));
                [~,ww] = contour(x_cont,y_cont,'LineColor','k');
                ww.LineWidth = 2;
            end
            title(sprintf('Component %i ',i),'fontsize',16,'fontweight','bold'); drawnow; %pause;
        else
            cla
            imagesc(reshape(b(:,i-nr),d1,d2)); axis equal; axis tight;
            title('Background component','fontsize',16,'fontweight','bold'); drawnow;
        end
        subplot(3,4,[10 11]);
        if i <= nr
            plot(1:T,DFF(i,:),'k','linewidth',1); hold all; plot(1:T,C(i,:),'b','linewidth',1);
            plot(1:T,3*S(i,:),'r','linewidth',1);
            if plot_df
                title(sprintf('Component %i (calcium DF/F value)',i),'fontsize',16,'fontweight','bold');
            else
                title(sprintf('Component %i (calcium raw value)',i),'fontsize',16,'fontweight','bold');
            end
%             leg = legend('Raw trace (filtered)','Inferred');
%             set(leg,'FontSize',14,'FontWeight','bold');
            drawnow;
            hold off;
        else
            plot(1:T,f(i-nr,:)); title('Background activity','fontsize',16,'fontweight','bold');
            drawnow;
            if make_gif
                frame = getframe(fig); %getframe(1);
                im = frame2im(frame);
                [imind,clm] = rgb2ind(im,256);
                if i == 1
                    imwrite(imind,clm,[name,'.gif'],'gif', 'Loopcount',inf);
                else
                    imwrite(imind,clm,[name,'.gif'],'gif','WriteMode','append');
                end

            end
        end
        
        
        
%         subplot(3,4,[3 4 7 8]);
%         set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');
%         if i <= nr
%            cla;
%            imagesc(mov(:,:,round(sld_t.Value))); hold on;
%            
%            [~,ww]=contour(x_cont, y_cont, 'LineColor','k');
%            ww.LineWidth=2;
%            title(sprintf('Frame %i ',round(sld_t.Value)),'fontsize',16,'fontweight','bold'); 
%            drawnow; %pause;
%            hold off;
%         end
        
        
            
    end

%     function rawDatUpdate(source,callbackdata)
%         t = source.Value;
%         plotRawDat(round(t));
%     end

%     function plotRawDat(t)
%         i = round(sld.Value);
%  
%         subplot(3,4,[3 4 7 8]);  
%         set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');
%         if i <= nr
%             cla;
%             imagesc('cData',mov(:,:,t),'XData',[0 796],'YData',[0 512]); hold on;
%             [~,ww]=contour(x_cont, y_cont, 'LineColor','k');
%             ww.LineWidth=2;
%             title(sprintf('Frame %i',t),'fontsize',16,'fontweight','bold');
%            
%            drawnow; %pause;
%            hold off;
%            
%         end 
%         
%       
%         
%     end

    function rejectComponent(src,event)
        r = round(sld.Value);
        display(sprintf('component %i rejected',r));
        rejectFile.reject(r,1) = 1;
    end

    function unrejectComponent(src,event)
        r = round(sld.Value);
        reject(r,1) = 0;
    end

    function KeyPressFcn(hObject,event,handles)
        
        h = get(gcf,'currentkey');
        switch h
%             case 'rightarrow'
% %                 disp('right')
%                 sld_t.Value = min(round(sld_t.Value)+1,T);
%                 plotRawDat(sld_t.Value);
%             case 'leftarrow'
% %                 disp('left')
%                 sld_t.Value = max(1, round(sld_t.Value-1));
%                 plotRawDat(sld_t.Value);
            case 'uparrow'
                sld.Value = min(round(sld.Value)+1,size(A,2));
                plot_component(sld.Value);
            case 'downarrow'
                sld.Value = max(1, round(sld.Value-1));
                plot_component(sld.Value);
                
            case 'r'
                reject(round(sld.Value)) = 1;
                sprintf('Rejected Component %i',round(sld.Value))
                
            case 'u'
                reject(sld.Value) = 0;
                sprintf('Unrejected Componenet %i', round(sld.Value))
                
        end
    end
            
end