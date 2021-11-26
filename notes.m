
% % Fjerner Bias fra MP (MP-Bias)    
% [m,n] = size(mean_table_MP1);
% real_MP1 = all_sat_MP1;
% real_MP2 = all_sat_MP2;
% for i = 1:m
%     real_MP1(mean_table_MP1(i,2):mean_table_MP1(i,3),mean_table_MP1(i,1)) = (real_MP1(mean_table_MP1(i,2):mean_table_MP1(i,3),mean_table_MP1(i,1))) - (mean_table_MP1(i,4));
%     real_MP2(mean_table_MP2(i,2):mean_table_MP2(i,3),mean_table_MP2(i,1)) = (real_MP2(mean_table_MP2(i,2):mean_table_MP2(i,3),mean_table_MP2(i,1))) - (mean_table_MP2(i,4));
% end
%     
% %Fjerner alle kolonner med kun "NaN" verdier til å plotte
% plot_MP1 = real_MP1(:,~all(isnan(real_MP1)));
% plot_MP2 = real_MP2(:,~all(isnan(real_MP2)));
% satelites = unique(mean_table_MP1(:,1));
% legend_list = append('PRN',string(satelites(:)));
% % Plot
% [m, n] = size(plot_MP1);
% [m, n] = size(plot_MP2);
% 
% % Test
% 
% MP1 = figure;
% for i = 1:n
%     plot(plot_MP1(:, i))
%     hold on
% end
% [row,kolonne] =  size(plot_MP1);
% mean_value = mean(abs(plot_MP1(1:row,kolonne)),'omitnan');
% tekst = [system, 'Estimated multipath on the P2-code(MP2)'];
% title(tekst);
% ylim([-10 10]);
% xlim([-5 7700]);
% legend(legend_list)
% xlabel('Epochs') 
% ylabel('Noise(meters)') 
% filnavn = append(system,'_MP1.png');
% exportgraphics(MP1,filnavn)
% MP2 = figure;
% for i = 1:n
%     plot(plot_MP2(:, i))
%     hold on
%     
% end
% 
% [row,kolonne] =  size(plot_MP2);
% mean_value = mean(abs(plot_MP2(1:row,kolonne)),'omitnan');
% tekst = [system, 'Estimated multipath on the P2-code(MP2)'];
% title(tekst);
% ylim([-10 10]);
% xlim([-5 7700]);
% xlabel('Epochs') 
% legend(legend_list)
% ylabel('Noise(meters)')
% filnavn = append(system,'_MP2.png');
% exportgraphics(MP2,filnavn)
% %plot(real_MP(:,4))

% [row,kolonne] =  size(mean_table_MP1);
% [nr_sat,m] = size(satelites)
% bias_table_MP1 = zeros(31,3);
% bias_table_MP2 = zeros(31,3);
% 
% for k = 1:row
%     bias_table_MP1(mean_table_MP1(k,1),:) = [mean_table_MP1(k,1),bias_table_MP1(mean_table_MP1(k,1),2)+ mean_table_MP1(k,3)-mean_table_MP1(k,2),bias_table_MP1(mean_table_MP1(k,1),3)+1];
%     bias_table_MP2(mean_table_MP2(k,1),:) = [mean_table_MP2(k,1),bias_table_MP2(mean_table_MP2(k,1),2)+ mean_table_MP2(k,3)-mean_table_MP2(k,2),bias_table_MP2(mean_table_MP2(k,1),3)+1];
% end
% bias_table_MP2 = bias_table_MP2(any(bias_table_MP2,2),:);
% bias_table_MP1 = bias_table_MP1(any(bias_table_MP1,2),:);
% 
% S_mp1 = [];
% S_mp2 = [];
% sum_mp = [];
% m_mp1 = []; % mean mp1 
% m_mp2 = []; % mean mp1 
% for k = 1:nr_sat
%     
%     rsum_mp1 = sum(abs(real_MP1(:,satelites(k))),'omitnan'); % sum mp1
%     rsum_mp2 = sum(abs(real_MP2(:,satelites(k))),'omitnan'); % sum mp2
%     sum_mp1 = sum(real_MP1(:,satelites(k)).^2,'omitnan');
%     sum_mp2 = sum(real_MP2(:,satelites(k)).^2,'omitnan');
%     nevner_mp1 = (bias_table_MP1(k,2)-bias_table_MP1(k,3));
%     nevner_mp2 = (bias_table_MP2(k,2)-bias_table_MP2(k,3));
%     m_mp1 = [m_mp1;satelites(k),rsum_mp1/nevner_mp1]; % mean mp1
%     m_mp2 = [m_mp2;satelites(k),rsum_mp2/nevner_mp2]; % mean mp2
%     S_mp1 = [S_mp1;satelites(k),sqrt(sum_mp1/nevner_mp1)];
%     S_mp2 = [S_mp2;satelites(k),sqrt(sum_mp2/(nevner_mp2))];
%end 



%          if ~any([phase1, phase2, code2] == 0)
%             Mp2 = code2 - (2*alfa/(alfa-1))*phase1 + (2*alfa/(alfa-1)-1)*phase2;                      
%             all_sat_MP2(epoch,SV) = Mp2;
% %          % get code observations (m)
% %          code1 = GNSS_obs{sysIndex}(SV, code1_index, epoch);
% %          code2 = GNSS_obs{sysIndex}(SV, code2_index, epoch);
% %          if ~any([phase1, phase2] == 0)
% %               % Fasebruddindikator
% %             IOD = (alfa/(alfa-1))*(phase1 - phase2); 
% %             fasebrudd(epoch, SV) = IOD;
% %             % Multipath + Bias
% %             Mp1 = code1 - (1 + 2/(alfa-1))*phase1 + (2/(alfa-1))*phase2;                   
% %             all_sat_MP1(epoch,SV) = Mp1;
% %             Mp2 = code2 - (2*alfa/(alfa-1))*phase1 + (2*alfa/(alfa-1)-1)*phase2;                      
% %             all_sat_MP2(epoch,SV) = Mp2;
%            
%          end   
%       end
%    end
% end