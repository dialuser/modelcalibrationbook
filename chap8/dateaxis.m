function dateaxis(tickaxis,dateform,startdate) 
%DATEAXIS Date axis labels. 
%	DATEAXIS(TICKAXIS,DATEFORM,STARTDATE) replaces axis tick labels with
%	date labels.  TICKAXIS determines which axis tick labels, X, Y, or Z,
%	should be converted.  The default TICKAXIS argument is 'x'.  DATEFORM
%	specifies which date format to use.  If no DATEFORM argument is
%	entered, this function determines the date format based on the span of
%	the axis limits. For example, if the difference between the axis
%	minimum and maximum  is less than 15, the tick labels will be converted
%	to 3 letter day of the week abbreviations.  STARTDATE determines which
%	date should be assigned to the first axis tick value.  The tick values
%	are treated as serial date numbers.  The default STARTDATE is the lower
%	axis limit converted to the appropriate date value.  For example, a
%	tick value of 1 is converted to the date 01-Jan-0000.  By entering
%	STARTDATE as '06-Apr-1995', the first tick value is assigned the date
%	April 6, 1995 and the axis tick labels will be set accordingly.
%      
%       DATEFORM    Format                Description 
% 
%          0        01-Mar-1995 15:45:17  (day-month-year, hour:minute) 
%          1        01-Mar-1995           (day-month-year) 
%          2        03/01/95              (month/day/year) 
%          3        Mar                   (month, three letter) 
%          4        M                     (month, single letter) 
%          5        3                     (month) 
%          6        03/01                 (month/day) 
%          7        1                     (day of month) 
%          8        Wed                   (day of week, three letter) 
%          9        W                     (day of week, single letter) 
%          10       1995                  (year, four digit) 
%          11       95                    (year, two digit) 
%          12       Mar95                 (month year) 
%          13       15:45:17              (hour:minute:second) 
%          14       03:45:17              (hour:minute:second AM or PM) 
%          15       15:45                 (hour:minute) 
%          16       03:45 PM              (hour:minute AM or PM) 
%          17       95/03/01              (year/month/day) 
% 
%	DATEAXIS('X') or DATEAXIS converts the X-axis labels to an 
%	automatically determined date format. 
% 
%	DATEAXIS('Y',6) converts the Y-axis labels to the 
%	month/day format. 
%     
%	DATEAXIS('X',2,'03/03/1995') converts the X-axis 
%	labels to the month/day/year format.  The minimum Xtick  
%	value is treated as March 3, 1995. 
 
%	Copyright 1995-2010 The MathWorks, Inc. 
%	$Revision: 1.11.2.4 $   $Date: 2010/10/08 16:42:35 $ 
 
dstr = []; 
 
switch nargin
  case 0 
    tickaxis = 'x'; 
    Lim = get(gca,[tickaxis,'lim']); 
    startdate = 0; 
  case 1 
    Lim = get(gca,[tickaxis,'lim']); 
    startdate = 0; 
  case 2   
    startdate = 0; 
  case 3 
    Lim = get(gca,[tickaxis,'lim']); 
    startdate = datenum(startdate)-Lim(1); 
end 
 
if nargin < 2 
  % Determine range of data and choose appropriate label format 
  Cond = Lim(2)-Lim(1); 
 
  if Cond <= 14 % Range less than 15 days, day of week   
    dateform = 7;  
  elseif Cond > 14 && Cond <= 31 % Range less than 32 days, day of month 
    dateform = 6; 
  elseif Cond > 31 && Cond <= 180 % Range less than 181 days, month/day 
    dateform = 5; 
  elseif Cond > 180 && Cond <= 365  % Range less than 366 days, 3 letter month 
    dateform = 3; 
  elseif Cond > 365 && Cond <= 365*3 % Range less than 3 years, month year  
    dateform = 11; 
  else % Range greater than 3 years, 2 digit year 
    dateform = 10; 
  end 
end 
 
% Get axis tick values and add appropriate start date. 
xl = get(gca,[tickaxis,'tick'])'+datenum(startdate); 
set(gca,[tickaxis,'tickmode'],'manual',[tickaxis,'limmode'],'manual') 
 
n = length(xl); 
 
% To guarantee that the day, month, and year strings have the 
% the same number of characters after the integer to string  
% conversion, add 100 to the day and month numbers and 10000 to 
% the year numbers for proper padding.  This also allows for 
% padding numbers with zeros.  For example, the first day of a  
% given month will be printed as 01 instead of 1.  The tick values 
% are converted into one long concatenated string by int2str. 
% Reshape matrix so each column is a single value plus the  
% appropriate padding number.  Transpose so each row is now this 
% single value and the needed columns can be extracted. 
 
switch dateform
  case {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16} 
    dstr = datestr(xl,dateform);
  case 17 % Year/Month/Day  (ISO format)    
    dstr = datestr(xl,25);
  otherwise
    error(message('finance:calendar:dateAxis'))       
end 
  
% Set axis tick labels 
set(gca,[tickaxis,'ticklabel'],dstr)
