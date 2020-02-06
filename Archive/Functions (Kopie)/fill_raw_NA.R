fill_raw_NA_time <- function(noe){
  #macht leere tibbles, die später befüllt werden können
  #Anzahl der Elemente der Zeitreihe = noe
  #Anzahl der Zeitreihen muss leider vorher bekannt sein
  data = tibble(time1 = numeric(noe),  time2 = numeric(noe),  time3 = numeric(noe),  time4 = numeric(noe),  time5 = numeric(noe),  time6 = numeric(noe),  time7 = numeric(noe),  time8 = numeric(noe),  time9 = numeric(noe),  time10 = numeric(noe),
                time11 = numeric(noe), time12 = numeric(noe), time13 = numeric(noe), time14 = numeric(noe), time15 = numeric(noe), time16 = numeric(noe), time17 = numeric(noe), time18 = numeric(noe), time19 = numeric(noe), time20 = numeric(noe),
                time21 = numeric(noe), time22 = numeric(noe), time23 = numeric(noe), time24 = numeric(noe), time25 = numeric(noe), time26 = numeric(noe), time27 = numeric(noe), time28 = numeric(noe), time29 = numeric(noe), time30 = numeric(noe), 
                time31 = numeric(noe), time32 = numeric(noe), time33 = numeric(noe), time34 = numeric(noe), time35 = numeric(noe), time36 = numeric(noe), time37 = numeric(noe), time38 = numeric(noe), time39 = numeric(noe), time40 = numeric(noe),
                time41 = numeric(noe), time42 = numeric(noe), time43 = numeric(noe), time44 = numeric(noe), time45 = numeric(noe), time46 = numeric(noe), time47 = numeric(noe), time48 = numeric(noe), time49 = numeric(noe), time50 = numeric(noe),
                time51 = numeric(noe), time52 = numeric(noe), time53 = numeric(noe), time54 = numeric(noe), time55 = numeric(noe), time56 = numeric(noe), time57 = numeric(noe), time58 = numeric(noe), time59 = numeric(noe), time60 = numeric(noe),
                time61 = numeric(noe), time62 = numeric(noe), time63 = numeric(noe), time64 = numeric(noe), time65 = numeric(noe), time66 = numeric(noe), time67 = numeric(noe), time68 = numeric(noe), time69 = numeric(noe), time70 = numeric(noe),
                time71 = numeric(noe), time72 = numeric(noe), time73 = numeric(noe), time74 = numeric(noe), time75 = numeric(noe), time76 = numeric(noe), time77 = numeric(noe), time78 = numeric(noe), time79 = numeric(noe), time80 = numeric(noe),
                time81 = numeric(noe), time82 = numeric(noe), time83 = numeric(noe), time84 = numeric(noe), time85 = numeric(noe), time86 = numeric(noe), time87 = numeric(noe), time88 = numeric(noe), time89 = numeric(noe), time90 = numeric(noe),
                time91 = numeric(noe), time92 = numeric(noe), time93 = numeric(noe), time94 = numeric(noe), time95 = numeric(noe), time96 = numeric(noe), time97 = numeric(noe), time98 = numeric(noe), time99 = numeric(noe), time100 = numeric(noe),
                time101 = numeric(noe), time102 = numeric(noe), time103 = numeric(noe), time104 = numeric(noe), time105 = numeric(noe), time106 = numeric(noe), time107 = numeric(noe), time108 = numeric(noe), time109 = numeric(noe), time110 = numeric(noe),
                time111 = numeric(noe), time112 = numeric(noe), time113 = numeric(noe), time114 = numeric(noe), time115 = numeric(noe), time116 = numeric(noe), time117 = numeric(noe), time118 = numeric(noe), time119 = numeric(noe), time120 = numeric(noe),
                time121 = numeric(noe), time122 = numeric(noe), time123 = numeric(noe), time124 = numeric(noe), time125 = numeric(noe), time126 = numeric(noe), time127 = numeric(noe), time128 = numeric(noe), time129 = numeric(noe), time130 = numeric(noe),
                time131 = numeric(noe), time132 = numeric(noe), time133 = numeric(noe), time134 = numeric(noe), time135 = numeric(noe), time136 = numeric(noe), time137 = numeric(noe), time138 = numeric(noe), time139 = numeric(noe), time140 = numeric(noe),
                time141 = numeric(noe), time142 = numeric(noe), time143 = numeric(noe), time144 = numeric(noe), time145 = numeric(noe), time146 = numeric(noe), time147 = numeric(noe), time148 = numeric(noe), time149 = numeric(noe), time150 = numeric(noe),
                time151 = numeric(noe), time152 = numeric(noe), time153 = numeric(noe), time154 = numeric(noe), time155 = numeric(noe), time156 = numeric(noe), time157 = numeric(noe), time158 = numeric(noe), time159 = numeric(noe), time160 = numeric(noe),
                time161 = numeric(noe), time162 = numeric(noe), time163 = numeric(noe), time164 = numeric(noe), time165 = numeric(noe), time166 = numeric(noe), time167 = numeric(noe), time168 = numeric(noe), time169 = numeric(noe), time170 = numeric(noe),
                time171 = numeric(noe), time172 = numeric(noe), time173 = numeric(noe), time174 = numeric(noe), time175 = numeric(noe), time176 = numeric(noe), time177 = numeric(noe), time178 = numeric(noe), time179 = numeric(noe), time180 = numeric(noe),
                time181 = numeric(noe), time182 = numeric(noe), time183 = numeric(noe), time184 = numeric(noe), time185 = numeric(noe), time186 = numeric(noe), time187 = numeric(noe), time188 = numeric(noe), time189 = numeric(noe), time190 = numeric(noe),
                time191 = numeric(noe), time192 = numeric(noe), time193 = numeric(noe), time194 = numeric(noe), time195 = numeric(noe), time196 = numeric(noe), time197 = numeric(noe), time198 = numeric(noe), time199 = numeric(noe), time200 = numeric(noe),
                time201 = numeric(noe), time202 = numeric(noe), time203 = numeric(noe), time204 = numeric(noe), time205 = numeric(noe), time206 = numeric(noe), time207 = numeric(noe), time208 = numeric(noe), time209 = numeric(noe), time210 = numeric(noe),
                time211 = numeric(noe), time212 = numeric(noe))
  return(data)
  }

fill_raw_NA_temp <- function(noe){
  #macht leere tibbles, die später befüllt werden können
  #Anzahl der Elemente der Zeitreihe = noe
  #Anzahl der Zeitreihen muss leider vorher bekannt sein
  data = tibble(temp1 = numeric(noe),  temp2 = numeric(noe),  temp3 = numeric(noe),  temp4 = numeric(noe),  temp5 = numeric(noe),  temp6 = numeric(noe),  temp7 = numeric(noe),  temp8 = numeric(noe),  temp9 = numeric(noe),  temp10 = numeric(noe),
                temp11 = numeric(noe), temp12 = numeric(noe), temp13 = numeric(noe), temp14 = numeric(noe), temp15 = numeric(noe), temp16 = numeric(noe), temp17 = numeric(noe), temp18 = numeric(noe), temp19 = numeric(noe), temp20 = numeric(noe),
                temp21 = numeric(noe), temp22 = numeric(noe), temp23 = numeric(noe), temp24 = numeric(noe), temp25 = numeric(noe), temp26 = numeric(noe), temp27 = numeric(noe), temp28 = numeric(noe), temp29 = numeric(noe), temp30 = numeric(noe), 
                temp31 = numeric(noe), temp32 = numeric(noe), temp33 = numeric(noe), temp34 = numeric(noe), temp35 = numeric(noe), temp36 = numeric(noe), temp37 = numeric(noe), temp38 = numeric(noe), temp39 = numeric(noe), temp40 = numeric(noe),
                temp41 = numeric(noe), temp42 = numeric(noe), temp43 = numeric(noe), temp44 = numeric(noe), temp45 = numeric(noe), temp46 = numeric(noe), temp47 = numeric(noe), temp48 = numeric(noe), temp49 = numeric(noe), temp50 = numeric(noe),
                temp51 = numeric(noe), temp52 = numeric(noe), temp53 = numeric(noe), temp54 = numeric(noe), temp55 = numeric(noe), temp56 = numeric(noe), temp57 = numeric(noe), temp58 = numeric(noe), temp59 = numeric(noe), temp60 = numeric(noe),
                temp61 = numeric(noe), temp62 = numeric(noe), temp63 = numeric(noe), temp64 = numeric(noe), temp65 = numeric(noe), temp66 = numeric(noe), temp67 = numeric(noe), temp68 = numeric(noe), temp69 = numeric(noe), temp70 = numeric(noe),
                temp71 = numeric(noe), temp72 = numeric(noe), temp73 = numeric(noe), temp74 = numeric(noe), temp75 = numeric(noe), temp76 = numeric(noe), temp77 = numeric(noe), temp78 = numeric(noe), temp79 = numeric(noe), temp80 = numeric(noe),
                temp81 = numeric(noe), temp82 = numeric(noe), temp83 = numeric(noe), temp84 = numeric(noe), temp85 = numeric(noe), temp86 = numeric(noe), temp87 = numeric(noe), temp88 = numeric(noe), temp89 = numeric(noe), temp90 = numeric(noe),
                temp91 = numeric(noe), temp92 = numeric(noe), temp93 = numeric(noe), temp94 = numeric(noe), temp95 = numeric(noe), temp96 = numeric(noe), temp97 = numeric(noe), temp98 = numeric(noe), temp99 = numeric(noe), temp100 = numeric(noe),
                temp101 = numeric(noe), temp102 = numeric(noe), temp103 = numeric(noe), temp104 = numeric(noe), temp105 = numeric(noe), temp106 = numeric(noe), temp107 = numeric(noe), temp108 = numeric(noe), temp109 = numeric(noe), temp110 = numeric(noe),
                temp111 = numeric(noe), temp112 = numeric(noe), temp113 = numeric(noe), temp114 = numeric(noe), temp115 = numeric(noe), temp116 = numeric(noe), temp117 = numeric(noe), temp118 = numeric(noe), temp119 = numeric(noe), temp120 = numeric(noe),
                temp121 = numeric(noe), temp122 = numeric(noe), temp123 = numeric(noe), temp124 = numeric(noe), temp125 = numeric(noe), temp126 = numeric(noe), temp127 = numeric(noe), temp128 = numeric(noe), temp129 = numeric(noe), temp130 = numeric(noe),
                temp131 = numeric(noe), temp132 = numeric(noe), temp133 = numeric(noe), temp134 = numeric(noe), temp135 = numeric(noe), temp136 = numeric(noe), temp137 = numeric(noe), temp138 = numeric(noe), temp139 = numeric(noe), temp140 = numeric(noe),
                temp141 = numeric(noe), temp142 = numeric(noe), temp143 = numeric(noe), temp144 = numeric(noe), temp145 = numeric(noe), temp146 = numeric(noe), temp147 = numeric(noe), temp148 = numeric(noe), temp149 = numeric(noe), temp150 = numeric(noe),
                temp151 = numeric(noe), temp152 = numeric(noe), temp153 = numeric(noe), temp154 = numeric(noe), temp155 = numeric(noe), temp156 = numeric(noe), temp157 = numeric(noe), temp158 = numeric(noe), temp159 = numeric(noe), temp160 = numeric(noe),
                temp161 = numeric(noe), temp162 = numeric(noe), temp163 = numeric(noe), temp164 = numeric(noe), temp165 = numeric(noe), temp166 = numeric(noe), temp167 = numeric(noe), temp168 = numeric(noe), temp169 = numeric(noe), temp170 = numeric(noe),
                temp171 = numeric(noe), temp172 = numeric(noe), temp173 = numeric(noe), temp174 = numeric(noe), temp175 = numeric(noe), temp176 = numeric(noe), temp177 = numeric(noe), temp178 = numeric(noe), temp179 = numeric(noe), temp180 = numeric(noe),
                temp181 = numeric(noe), temp182 = numeric(noe), temp183 = numeric(noe), temp184 = numeric(noe), temp185 = numeric(noe), temp186 = numeric(noe), temp187 = numeric(noe), temp188 = numeric(noe), temp189 = numeric(noe), temp190 = numeric(noe),
                temp191 = numeric(noe), temp192 = numeric(noe), temp193 = numeric(noe), temp194 = numeric(noe), temp195 = numeric(noe), temp196 = numeric(noe), temp197 = numeric(noe), temp198 = numeric(noe), temp199 = numeric(noe), temp200 = numeric(noe),
                temp201 = numeric(noe), temp202 = numeric(noe), temp203 = numeric(noe), temp204 = numeric(noe), temp205 = numeric(noe), temp206 = numeric(noe), temp207 = numeric(noe), temp208 = numeric(noe), temp209 = numeric(noe), temp210 = numeric(noe),
                temp211 = numeric(noe), temp212 = numeric(noe))
  return(data)
}

fill_raw_NA_prec <- function(noe){
  #macht leere tibbles, die später befüllt werden können
  #Anzahl der Elemente der Zeitreihe = noe
  #Anzahl der Zeitreihen muss leider vorher bekannt sein
  data = tibble(prec1 = numeric(noe),  prec2 = numeric(noe),  prec3 = numeric(noe),  prec4 = numeric(noe),  prec5 = numeric(noe),  prec6 = numeric(noe),  prec7 = numeric(noe),  prec8 = numeric(noe),  prec9 = numeric(noe),  prec10 = numeric(noe),
                prec11 = numeric(noe), prec12 = numeric(noe), prec13 = numeric(noe), prec14 = numeric(noe), prec15 = numeric(noe), prec16 = numeric(noe), prec17 = numeric(noe), prec18 = numeric(noe), prec19 = numeric(noe), prec20 = numeric(noe),
                prec21 = numeric(noe), prec22 = numeric(noe), prec23 = numeric(noe), prec24 = numeric(noe), prec25 = numeric(noe), prec26 = numeric(noe), prec27 = numeric(noe), prec28 = numeric(noe), prec29 = numeric(noe), prec30 = numeric(noe), 
                prec31 = numeric(noe), prec32 = numeric(noe), prec33 = numeric(noe), prec34 = numeric(noe), prec35 = numeric(noe), prec36 = numeric(noe), prec37 = numeric(noe), prec38 = numeric(noe), prec39 = numeric(noe), prec40 = numeric(noe),
                prec41 = numeric(noe), prec42 = numeric(noe), prec43 = numeric(noe), prec44 = numeric(noe), prec45 = numeric(noe), prec46 = numeric(noe), prec47 = numeric(noe), prec48 = numeric(noe), prec49 = numeric(noe), prec50 = numeric(noe),
                prec51 = numeric(noe), prec52 = numeric(noe), prec53 = numeric(noe), prec54 = numeric(noe), prec55 = numeric(noe), prec56 = numeric(noe), prec57 = numeric(noe), prec58 = numeric(noe), prec59 = numeric(noe), prec60 = numeric(noe),
                prec61 = numeric(noe), prec62 = numeric(noe), prec63 = numeric(noe), prec64 = numeric(noe), prec65 = numeric(noe), prec66 = numeric(noe), prec67 = numeric(noe), prec68 = numeric(noe), prec69 = numeric(noe), prec70 = numeric(noe),
                prec71 = numeric(noe), prec72 = numeric(noe), prec73 = numeric(noe), prec74 = numeric(noe), prec75 = numeric(noe), prec76 = numeric(noe), prec77 = numeric(noe), prec78 = numeric(noe), prec79 = numeric(noe), prec80 = numeric(noe),
                prec81 = numeric(noe), prec82 = numeric(noe), prec83 = numeric(noe), prec84 = numeric(noe), prec85 = numeric(noe), prec86 = numeric(noe), prec87 = numeric(noe), prec88 = numeric(noe), prec89 = numeric(noe), prec90 = numeric(noe),
                prec91 = numeric(noe), prec92 = numeric(noe), prec93 = numeric(noe), prec94 = numeric(noe), prec95 = numeric(noe), prec96 = numeric(noe), prec97 = numeric(noe), prec98 = numeric(noe), prec99 = numeric(noe), prec100 = numeric(noe),
                prec101 = numeric(noe), prec102 = numeric(noe), prec103 = numeric(noe), prec104 = numeric(noe), prec105 = numeric(noe), prec106 = numeric(noe), prec107 = numeric(noe), prec108 = numeric(noe), prec109 = numeric(noe), prec110 = numeric(noe),
                prec111 = numeric(noe), prec112 = numeric(noe), prec113 = numeric(noe), prec114 = numeric(noe), prec115 = numeric(noe), prec116 = numeric(noe), prec117 = numeric(noe), prec118 = numeric(noe), prec119 = numeric(noe), prec120 = numeric(noe),
                prec121 = numeric(noe), prec122 = numeric(noe), prec123 = numeric(noe), prec124 = numeric(noe), prec125 = numeric(noe), prec126 = numeric(noe), prec127 = numeric(noe), prec128 = numeric(noe), prec129 = numeric(noe), prec130 = numeric(noe),
                prec131 = numeric(noe), prec132 = numeric(noe), prec133 = numeric(noe), prec134 = numeric(noe), prec135 = numeric(noe), prec136 = numeric(noe), prec137 = numeric(noe), prec138 = numeric(noe), prec139 = numeric(noe), prec140 = numeric(noe),
                prec141 = numeric(noe), prec142 = numeric(noe), prec143 = numeric(noe), prec144 = numeric(noe), prec145 = numeric(noe), prec146 = numeric(noe), prec147 = numeric(noe), prec148 = numeric(noe), prec149 = numeric(noe), prec150 = numeric(noe),
                prec151 = numeric(noe), prec152 = numeric(noe), prec153 = numeric(noe), prec154 = numeric(noe), prec155 = numeric(noe), prec156 = numeric(noe), prec157 = numeric(noe), prec158 = numeric(noe), prec159 = numeric(noe), prec160 = numeric(noe),
                prec161 = numeric(noe), prec162 = numeric(noe), prec163 = numeric(noe), prec164 = numeric(noe), prec165 = numeric(noe), prec166 = numeric(noe), prec167 = numeric(noe), prec168 = numeric(noe), prec169 = numeric(noe), prec170 = numeric(noe),
                prec171 = numeric(noe), prec172 = numeric(noe), prec173 = numeric(noe), prec174 = numeric(noe), prec175 = numeric(noe), prec176 = numeric(noe), prec177 = numeric(noe), prec178 = numeric(noe), prec179 = numeric(noe), prec180 = numeric(noe),
                prec181 = numeric(noe), prec182 = numeric(noe), prec183 = numeric(noe), prec184 = numeric(noe), prec185 = numeric(noe), prec186 = numeric(noe), prec187 = numeric(noe), prec188 = numeric(noe), prec189 = numeric(noe), prec190 = numeric(noe),
                prec191 = numeric(noe), prec192 = numeric(noe), prec193 = numeric(noe), prec194 = numeric(noe), prec195 = numeric(noe), prec196 = numeric(noe), prec197 = numeric(noe), prec198 = numeric(noe), prec199 = numeric(noe), prec200 = numeric(noe),
                prec201 = numeric(noe), prec202 = numeric(noe), prec203 = numeric(noe), prec204 = numeric(noe), prec205 = numeric(noe), prec206 = numeric(noe), prec207 = numeric(noe), prec208 = numeric(noe), prec209 = numeric(noe), prec210 = numeric(noe),
                prec211 = numeric(noe), prec212 = numeric(noe))
  return(data)
}

fill_raw_NA_isot <- function(noe){
  #macht leere tibbles, die später befüllt werden können
  #Anzahl der Elemente der Zeitreihe = noe
  #Anzahl der Zeitreihen muss leider vorher bekannt sein
  data = tibble(isot1 = numeric(noe),  isot2 = numeric(noe),  isot3 = numeric(noe),  isot4 = numeric(noe),  isot5 = numeric(noe),  isot6 = numeric(noe),  isot7 = numeric(noe),  isot8 = numeric(noe),  isot9 = numeric(noe),  isot10 = numeric(noe),
                isot11 = numeric(noe), isot12 = numeric(noe), isot13 = numeric(noe), isot14 = numeric(noe), isot15 = numeric(noe), isot16 = numeric(noe), isot17 = numeric(noe), isot18 = numeric(noe), isot19 = numeric(noe), isot20 = numeric(noe),
                isot21 = numeric(noe), isot22 = numeric(noe), isot23 = numeric(noe), isot24 = numeric(noe), isot25 = numeric(noe), isot26 = numeric(noe), isot27 = numeric(noe), isot28 = numeric(noe), isot29 = numeric(noe), isot30 = numeric(noe), 
                isot31 = numeric(noe), isot32 = numeric(noe), isot33 = numeric(noe), isot34 = numeric(noe), isot35 = numeric(noe), isot36 = numeric(noe), isot37 = numeric(noe), isot38 = numeric(noe), isot39 = numeric(noe), isot40 = numeric(noe),
                isot41 = numeric(noe), isot42 = numeric(noe), isot43 = numeric(noe), isot44 = numeric(noe), isot45 = numeric(noe), isot46 = numeric(noe), isot47 = numeric(noe), isot48 = numeric(noe), isot49 = numeric(noe), isot50 = numeric(noe),
                isot51 = numeric(noe), isot52 = numeric(noe), isot53 = numeric(noe), isot54 = numeric(noe), isot55 = numeric(noe), isot56 = numeric(noe), isot57 = numeric(noe), isot58 = numeric(noe), isot59 = numeric(noe), isot60 = numeric(noe),
                isot61 = numeric(noe), isot62 = numeric(noe), isot63 = numeric(noe), isot64 = numeric(noe), isot65 = numeric(noe), isot66 = numeric(noe), isot67 = numeric(noe), isot68 = numeric(noe), isot69 = numeric(noe), isot70 = numeric(noe),
                isot71 = numeric(noe), isot72 = numeric(noe), isot73 = numeric(noe), isot74 = numeric(noe), isot75 = numeric(noe), isot76 = numeric(noe), isot77 = numeric(noe), isot78 = numeric(noe), isot79 = numeric(noe), isot80 = numeric(noe),
                isot81 = numeric(noe), isot82 = numeric(noe), isot83 = numeric(noe), isot84 = numeric(noe), isot85 = numeric(noe), isot86 = numeric(noe), isot87 = numeric(noe), isot88 = numeric(noe), isot89 = numeric(noe), isot90 = numeric(noe),
                isot91 = numeric(noe), isot92 = numeric(noe), isot93 = numeric(noe), isot94 = numeric(noe), isot95 = numeric(noe), isot96 = numeric(noe), isot97 = numeric(noe), isot98 = numeric(noe), isot99 = numeric(noe), isot100 = numeric(noe),
                isot101 = numeric(noe), isot102 = numeric(noe), isot103 = numeric(noe), isot104 = numeric(noe), isot105 = numeric(noe), isot106 = numeric(noe), isot107 = numeric(noe), isot108 = numeric(noe), isot109 = numeric(noe), isot110 = numeric(noe),
                isot111 = numeric(noe), isot112 = numeric(noe), isot113 = numeric(noe), isot114 = numeric(noe), isot115 = numeric(noe), isot116 = numeric(noe), isot117 = numeric(noe), isot118 = numeric(noe), isot119 = numeric(noe), isot120 = numeric(noe),
                isot121 = numeric(noe), isot122 = numeric(noe), isot123 = numeric(noe), isot124 = numeric(noe), isot125 = numeric(noe), isot126 = numeric(noe), isot127 = numeric(noe), isot128 = numeric(noe), isot129 = numeric(noe), isot130 = numeric(noe),
                isot131 = numeric(noe), isot132 = numeric(noe), isot133 = numeric(noe), isot134 = numeric(noe), isot135 = numeric(noe), isot136 = numeric(noe), isot137 = numeric(noe), isot138 = numeric(noe), isot139 = numeric(noe), isot140 = numeric(noe),
                isot141 = numeric(noe), isot142 = numeric(noe), isot143 = numeric(noe), isot144 = numeric(noe), isot145 = numeric(noe), isot146 = numeric(noe), isot147 = numeric(noe), isot148 = numeric(noe), isot149 = numeric(noe), isot150 = numeric(noe),
                isot151 = numeric(noe), isot152 = numeric(noe), isot153 = numeric(noe), isot154 = numeric(noe), isot155 = numeric(noe), isot156 = numeric(noe), isot157 = numeric(noe), isot158 = numeric(noe), isot159 = numeric(noe), isot160 = numeric(noe),
                isot161 = numeric(noe), isot162 = numeric(noe), isot163 = numeric(noe), isot164 = numeric(noe), isot165 = numeric(noe), isot166 = numeric(noe), isot167 = numeric(noe), isot168 = numeric(noe), isot169 = numeric(noe), isot170 = numeric(noe),
                isot171 = numeric(noe), isot172 = numeric(noe), isot173 = numeric(noe), isot174 = numeric(noe))
  return(data)
}
