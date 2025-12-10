#ifndef SRC_CORE_INCLUDE_MATH_Z_CONSTANTS_H_
#define SRC_CORE_INCLUDE_MATH_Z_CONSTANTS_H_

#include "math/hal/bigfixedpoint.h"

namespace lbcrypto {

// For X^32 - X + 2
const std::vector<BigComplex> z_upper_roots_32 = {
    BigComplex(BigFixedPoint(BigInteger("350551088990192644005587638168552186653"), 128, true),
               BigFixedPoint(BigInteger("34900550222895179714397787198203513173"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("336695400736306456376811924693546769860"), 128, true),
               BigFixedPoint(BigInteger("103299331694394466680808556291558998110"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("309551080397728457701091292115470551977"), 128, true),
               BigFixedPoint(BigInteger("167549268673586454371916803328721791984"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("270228898653084472294462899074132324165"), 128, true),
               BigFixedPoint(BigInteger("225075857771326758427915586988199792655"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("220337572168790987178561844883347289138"), 128, true),
               BigFixedPoint(BigInteger("273583686243471418078046036292938539652"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("161917343781254382163963077301550563567"), 128, true),
               BigFixedPoint(BigInteger("311152443074692736824352186014335408156"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("97355433476106157288077219655957918362"), 128, true),
               BigFixedPoint(BigInteger("336318161343136424694695901841775800828"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("29286437520563535328487200882786010496"), 128, true),
               BigFixedPoint(BigInteger("348136658982339264967867138036106806354"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("39518997925195344591398501225959432883"), 128, false),
               BigFixedPoint(BigInteger("346226947218601574485487145087988588362"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("106274150951931075782750588501284147648"), 128, false),
               BigFixedPoint(BigInteger("330793168550893306693876614221731825723"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("168302822651450193478236284374557061220"), 128, false),
               BigFixedPoint(BigInteger("302624017298379579243622911980893710731"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("223172712788269601050257487780680679410"), 128, false),
               BigFixedPoint(BigInteger("263067441411603613475403140531675225854"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("268838277499134743505006047266580888328"), 128, false),
               BigFixedPoint(BigInteger("213972753153500930100328268852044496823"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("303792157693636518457724129356882969974"), 128, false),
               BigFixedPoint(BigInteger("157576849671356512677734982596091307782"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("327191119747388911900252417088460946142"), 128, false),
               BigFixedPoint(BigInteger("96289389574348118233366899415301954126"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("338833016467020703571417641180937488614"), 128, false),
               BigFixedPoint(BigInteger("32362949412358133459608716350368962602"), 128, false)),
};

const auto z_upper_roots_32_scale = 128;

const auto z_upper_roots       = z_upper_roots_32;
const auto z_upper_roots_scale = z_upper_roots_32_scale;
const size_t zN                = 32;

BigCMatrix getZU();

BigCMatrix getZUInverse();

const BigCVector Z_ROOTS_N8 = {
    BigComplex(BigFixedPoint(BigInteger("359124072811535186610690853057402903661"), 128, true),
               BigFixedPoint(BigInteger("156682980275126563459473712549117854603"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("127942916012362888167921749339779666051"), 128, true),
               BigFixedPoint(BigInteger("361530722988969544481784436970553558351"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("163650965418420061481525854130989829121"), 128, false),
               BigFixedPoint(BigInteger("327582167183704011019777000705434796086"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("323416023405478013297086748266192740590"), 128, false),
               BigFixedPoint(BigInteger("118974601821943721998364824569673250474"), 128, false)),
};
const BigCVector Z_ROOTS_N16 = {
    BigComplex(BigFixedPoint(BigInteger("357580614754354758148531447970617000661"), 128, true),
               BigFixedPoint(BigInteger("72754167196966292104314199782731506654"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("300047663324999516624210886330771269624"), 128, true),
               BigFixedPoint(BigInteger("206155599878541348717138407546612184696"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("194911227399247327269977016223988976357"), 128, true),
               BigFixedPoint(BigInteger("305382816964920193167762588519087191920"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("60282570591716603101266424038243067083"), 128, true),
               BigFixedPoint(BigInteger("354456894667915792617878701713209385694"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("80784262001292603479136660970669486879"), 128, false),
               BigFixedPoint(BigInteger("346445020591532988638702568080557856406"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("204637768854078910624306365031060162530"), 128, false),
               BigFixedPoint(BigInteger("284972432336483361176802995403837422402"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("292276954111974733171082649947546090607"), 128, false),
               BigFixedPoint(BigInteger("183819236581063928453505381579739081840"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("335123091102971957869460098614344573710"), 128, false),
               BigFixedPoint(BigInteger("62755734557348253032878164961445143012"), 128, false)),
};
const BigCVector Z_ROOTS_N32 = {
    BigComplex(BigFixedPoint(BigInteger("350551088990192644005587638168552186653"), 128, true),
               BigFixedPoint(BigInteger("34900550222895179714397787198203513173"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("336695400736306456376811924693546769860"), 128, true),
               BigFixedPoint(BigInteger("103299331694394466680808556291558998110"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("309551080397728457701091292115470551977"), 128, true),
               BigFixedPoint(BigInteger("167549268673586454371916803328721791984"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("270228898653084472294462899074132324165"), 128, true),
               BigFixedPoint(BigInteger("225075857771326758427915586988199792655"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("220337572168790987178561844883347289138"), 128, true),
               BigFixedPoint(BigInteger("273583686243471418078046036292938539652"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("161917343781254382163963077301550563567"), 128, true),
               BigFixedPoint(BigInteger("311152443074692736824352186014335408156"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("97355433476106157288077219655957918362"), 128, true),
               BigFixedPoint(BigInteger("336318161343136424694695901841775800828"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("29286437520563535328487200882786010496"), 128, true),
               BigFixedPoint(BigInteger("348136658982339264967867138036106806354"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("39518997925195344591398501225959432883"), 128, false),
               BigFixedPoint(BigInteger("346226947218601574485487145087988588362"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("106274150951931075782750588501284147648"), 128, false),
               BigFixedPoint(BigInteger("330793168550893306693876614221731825723"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("168302822651450193478236284374557061220"), 128, false),
               BigFixedPoint(BigInteger("302624017298379579243622911980893710731"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("223172712788269601050257487780680679410"), 128, false),
               BigFixedPoint(BigInteger("263067441411603613475403140531675225854"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("268838277499134743505006047266580888328"), 128, false),
               BigFixedPoint(BigInteger("213972753153500930100328268852044496823"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("303792157693636518457724129356882969974"), 128, false),
               BigFixedPoint(BigInteger("157576849671356512677734982596091307782"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("327191119747388911900252417088460946142"), 128, false),
               BigFixedPoint(BigInteger("96289389574348118233366899415301954126"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("338833016467020703571417641180937488614"), 128, false),
               BigFixedPoint(BigInteger("32362949412358133459608716350368962602"), 128, false)),
};
const BigCVector Z_ROOTS_N64 = {
    BigComplex(BigFixedPoint(BigInteger("345782382768622893681963401448427898029"), 128, true),
               BigFixedPoint(BigInteger("17077306551372370810130247289563054112"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("342405286167006616104272679186613693604"), 128, true),
               BigFixedPoint(BigInteger("51063768767812853759267517276440925613"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("335684653456421547424442279525308350202"), 128, true),
               BigFixedPoint(BigInteger("84547459685791750395322149498754820193"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("325687271314630161660208052480672211255"), 128, true),
               BigFixedPoint(BigInteger("117198786494194762481512155704845056237"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("312512487251139019726072948662562208179"), 128, true),
               BigFixedPoint(BigInteger("148696480955218424558925142747080310379"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("296291219693069351048069937932644278857"), 128, true),
               BigFixedPoint(BigInteger("178730814836784585934091041688973274301"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("277184653141781729064681383692073336510"), 128, true),
               BigFixedPoint(BigInteger("207006702389922155776633147023636044372"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("255382631061828761393939614760687633559"), 128, true),
               BigFixedPoint(BigInteger("233246659645471238489308812049612276822"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("231101762036653735641835465279513483112"), 128, true),
               BigFixedPoint(BigInteger("257193591817794269088334104944013154924"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("204583257397787532099692355691355637319"), 128, true),
               BigFixedPoint(BigInteger("278613381928043479536542640077085181602"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("176090520958796792151684077165811983531"), 128, true),
               BigFixedPoint(BigInteger("297297255880365950935670182123312539677"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("145906513608118335885148193446290158527"), 128, true),
               BigFixedPoint(BigInteger("313063901621348585261744602459211782193"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("114330917272280526697576002156849896742"), 128, true),
               BigFixedPoint(BigInteger("325761322659738478048112596752588245545"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("81677124073488505954567391490098805922"), 128, true),
               BigFixedPoint(BigInteger("335268409083212843063340500613138477897"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("48269077270658794175472468304454533424"), 128, true),
               BigFixedPoint(BigInteger("341496212226481724524825912030575708205"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("14437990654826877108730547908067699438"), 128, true),
               BigFixedPoint(BigInteger("344388912231851359245712019344600207568"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("19481027713614270916657575720640649881"), 128, false),
               BigFixedPoint(BigInteger("343924470750769753248434163202513630816"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("53152423427826865766671048674191534627"), 128, false),
               BigFixedPoint(BigInteger("340114963704569587556699087847313256291"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("86243694617983255476829546911083929198"), 128, false),
               BigFixedPoint(BigInteger("333006590899766137860699657198694469622"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("118428926139984546657028007091204146214"), 128, false),
               BigFixedPoint(BigInteger("322679359576025494049286154563749908809"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("149392350458552483266058810791615997253"), 128, false),
               BigFixedPoint(BigInteger("309246436256479502314147527923531273152"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("178831935497824401290979186243808163315"), 128, false),
               BigFixedPoint(BigInteger("292853153202281558238070494053709348994"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("206463009072511083318194837251997731372"), 128, false),
               BigFixedPoint(BigInteger("273675638529907370012373951271030073982"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("232021934933830847250170168915294584786"), 128, false),
               BigFixedPoint(BigInteger("251919007047156328029803905526024276684"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("255269846174124869224583281394161026572"), 128, false),
               BigFixedPoint(BigInteger("227814995517252565325099776446496832246"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("275996395032360381148372509665527184202"), 128, false),
               BigFixedPoint(BigInteger("201618848130921411989812718440959519845"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("294023353999451160759043762351241867415"), 128, false),
               BigFixedPoint(BigInteger("173605168679034640797811756741931171060"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("309207649510494193856307686851957832063"), 128, false),
               BigFixedPoint(BigInteger("144062413573597809304023039687362531072"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("321443006846364342033233440671515137220"), 128, false),
               BigFixedPoint(BigInteger("113285847977699889107009096636872826161"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("330658968359936235109335002810174947792"), 128, false),
               BigFixedPoint(BigInteger("81569354560965542887680311904837222556"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("336816073143344975132084641642573087165"), 128, false),
               BigFixedPoint(BigInteger("49197604931280020061491849020858663676"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("339897153198907268612807292144443989133"), 128, false),
               BigFixedPoint(BigInteger("16441351304266807328870060951377273088"), 128, false)),
};
const BigCVector Z_ROOTS_N128 = {
    BigComplex(BigFixedPoint(BigInteger("343119165158351502202114921083918160607"), 128, true),
               BigFixedPoint(BigInteger("8445230262408170414710525911923979863"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("342286768299261249663092570720758204394"), 128, true),
               BigFixedPoint(BigInteger("25315119529903911884458726357125114939"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("340624011557953935681717439314745241930"), 128, true),
               BigFixedPoint(BigInteger("42123345521038082466654701457161616909"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("338134963897962496855538580572822899388"), 128, true),
               BigFixedPoint(BigInteger("58828967580866093795453906946573148931"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("334825716301157647535838535286922026820"), 128, true),
               BigFixedPoint(BigInteger("75391296919231694003989896991759267197"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("330704366845495045037151291913277532226"), 128, true),
               BigFixedPoint(BigInteger("91769996499480233872926784948326276727"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("325781000865017592270329428252011717270"), 128, true),
               BigFixedPoint(BigInteger("107925180072558935257495858760117596435"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("320067666240210658487505452969952890233"), 128, true),
               BigFixedPoint(BigInteger("123817510116021497411868184990925151288"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("313578343878591003450581379611289917832"), 128, true),
               BigFixedPoint(BigInteger("139408294440212911354764247324732797830"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("306328913457032201047618987833879966890"), 128, true),
               BigFixedPoint(BigInteger("154659581227269900264896765460915621431"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("298337114508758100400453600013616292626"), 128, true),
               BigFixedPoint(BigInteger("169534252272525810144260849910394470311"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("289622502949138963597583606781097337917"), 128, true),
               BigFixedPoint(BigInteger("183996114202447615972758498524229370252"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("280206403145369843732166807079314777349"), 128, true),
               BigFixedPoint(BigInteger("198009987448347300510880081796471256097"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("270111855645764723659388423786344916906"), 128, true),
               BigFixedPoint(BigInteger("211541792760789265063476337514300772622"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("259363560694729752913048412460909024151"), 128, true),
               BigFixedPoint(BigInteger("224558635055847439872375641289318967820"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("247987817669450879719683635884186985307"), 128, true),
               BigFixedPoint(BigInteger("237028884391136887515224946447466263972"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("236012460583910798620649287317304915401"), 128, true),
               BigFixedPoint(BigInteger("248922253876840117758636440944953065511"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("223466789815001908040388964706695414132"), 128, true),
               BigFixedPoint(BigInteger("260209874334751832022726256928096971727"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("210381500214188999129279261079302557594"), 128, true),
               BigFixedPoint(BigInteger("270864365526659633484305196559487150467"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("196788605776358960960103257311199627319"), 128, true),
               BigFixedPoint(BigInteger("280859903782142894871270335935082348731"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("182721361045133810044194141105554388121"), 128, true),
               BigFixedPoint(BigInteger("290172285865085981907510172402458343117"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("168214179440973727156191620091633926513"), 128, true),
               BigFixedPoint(BigInteger("298778988927841476340336626989804273643"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("153302548704810557967941831822086848154"), 128, true),
               BigFixedPoint(BigInteger("306659226412017007130942325233003065835"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("138022943655676614498723518228135261533"), 128, true),
               BigFixedPoint(BigInteger("313793999765265016653713358455432136714"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("122412736465769782021312496466204460593"), 128, true),
               BigFixedPoint(BigInteger("320166145854192489432691849969736624505"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("106510104660557629344640838719047424881"), 128, true),
               BigFixedPoint(BigInteger("325760379964534964613719657084857401518"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("90353937054795083638971543748599743106"), 128, true),
               BigFixedPoint(BigInteger("330563334291004810026653886940379494936"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("73983737837625920312150545433941403406"), 128, true),
               BigFixedPoint(BigInteger("334563591830664636372771428859118386609"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("57439529021158392265046202884079625678"), 128, true),
               BigFixedPoint(BigInteger("337751715605213703653159176391608565962"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("40761751466934991263283827849819053325"), 128, true),
               BigFixedPoint(BigInteger("340120273149107332605496053794316435023"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("23991164703423224491298212885922798015"), 128, true),
               BigFixedPoint(BigInteger("341663856211826410662201859834252049854"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("7168745744886518906170578802686823632"), 128, true),
               BigFixedPoint(BigInteger("342379095633706073698399579805036011555"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("9664412882420572195505503515950721605"), 128, false),
               BigFixedPoint(BigInteger("342264671365295937825809381654398380388"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("26467205707043382480663317811085399727"), 128, false),
               BigFixedPoint(BigInteger("341321317609963104344051463400047113356"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("43198617280718459902127166603572442038"), 128, false),
               BigFixedPoint(BigInteger("339551823077971168569647656552525157147"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("59817825367548669930797896168937062920"), 128, false),
               BigFixedPoint(BigInteger("336961026347052423954971845763192998630"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("76284304186188196738703739972266548978"), 128, false),
               BigFixedPoint(BigInteger("333555806328844710265349847818113035160"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("92557927552946780147811086415088996988"), 128, false),
               BigFixedPoint(BigInteger("329345067841573792421017057388546480745"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("108599071792121324996846648779935561547"), 128, false),
               BigFixedPoint(BigInteger("324339722285822322058770394641619880101"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("124368718300875824341570692574485132833"), 128, false),
               BigFixedPoint(BigInteger("318552663410564500140328365702345735109"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("139828555679206386726577856118867829312"), 128, false),
               BigFixedPoint(BigInteger("311998738138830585692557668366468987826"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("154941081360023128304249590287033597456"), 128, false),
               BigFixedPoint(BigInteger("304694712393815816106822881020357166601"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("169669702698413503053010414379533075987"), 128, false),
               BigFixedPoint(BigInteger("296659231823754633658924800289078700216"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("183978837499720447393345494708391871913"), 128, false),
               BigFixedPoint(BigInteger("287912777263579372424814413697594190415"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("197834013978285427526102510972766468370"), 128, false),
               BigFixedPoint(BigInteger("278477614688842859073841504912527576142"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("211201970134947684690158725239537370877"), 128, false),
               BigFixedPoint(BigInteger("268377739307924628490103583229265117638"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("224050752510206950670585443966104240133"), 128, false),
               BigFixedPoint(BigInteger("257638813297942570128996025857755621263"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("236349814194783876855068405485811974799"), 128, false),
               BigFixedPoint(BigInteger("246288096515689630834328200445362532134"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("248070111837172249812387417016338888752"), 128, false),
               BigFixedPoint(BigInteger("234354369309232219917440033982268642828"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("259184201148413721231376575389733491890"), 128, false),
               BigFixedPoint(BigInteger("221867846328717373013987696722643718330"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("269666330030679809711325331406777953628"), 128, false),
               BigFixedPoint(BigInteger("208860080010823140219396833753313499133"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("279492527907850282335386765327283334785"), 128, false),
               BigFixedPoint(BigInteger("195363852236917147277964397479733215584"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("288640689078021859528658474008340735368"), 128, false),
               BigFixedPoint(BigInteger("181413052619307322413944328840368626262"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("297090646927351870385752316892612352733"), 128, false),
               BigFixedPoint(BigInteger("167042542073162554808553030532077621548"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("304824234681108455232318090394303103378"), 128, false),
               BigFixedPoint(BigInteger("152288000948007407741949539443009105867"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("311825327153421537410378963364161599985"), 128, false),
               BigFixedPoint(BigInteger("137185762217133881831875407669521642303"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("318079856963936378321882495661345975392"), 128, false),
               BigFixedPoint(BigInteger("121772632237687761284097786556524043190"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("323575798361587282791140466505034657315"), 128, false),
               BigFixedPoint(BigInteger("106085704478537868474216152081989122152"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("328303112729329796099380426878087961407"), 128, false),
               BigFixedPoint(BigInteger("90162175208786663650238985669317406057"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("332253652662570190770048483646863738942"), 128, false),
               BigFixedPoint(BigInteger("74039173894483011890281313248589892199"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("335421026601794345596031424494531343321"), 128, false),
               BigFixedPoint(BigInteger("57753623925607795577494784534213911235"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("337800433118522047782969176243406131288"), 128, false),
               BigFixedPoint(BigInteger("41342149864611275157683690531707186982"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("339388481901658043426394317558112104426"), 128, false),
               BigFixedPoint(BigInteger("24841044275025775210360173562719093220"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("340183025076584028525603984230960495155"), 128, false),
               BigFixedPoint(BigInteger("8286299697137357239525899556631059252"), 128, false)),
};

const std::map<uint32_t, BigCVector> Z_ROOTS_MAP = {
    {8, Z_ROOTS_N8}, {16, Z_ROOTS_N16}, {32, Z_ROOTS_N32}, {64, Z_ROOTS_N64}, {128, Z_ROOTS_N128},
};

BigCMatrix GetZU(uint32_t zN);

BigCMatrix GetZUInverse(uint32_t zN);

// Now find the primitive roots for x^N+1

// this is zeta_0 to zeta_15 for N=32
// other parts are just conjugates
// we use the zeta notataion as in original CKKS paper
const std::vector<BigComplex> r_roots_32 = {
    BigComplex(BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, false),
               BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("300102255270364911803307903244229695504"), 128, false),
               BigFixedPoint(BigInteger("160407997365957196655387116530659854223"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("263041826724899848421390153303749767787"), 128, true),
               BigFixedPoint(BigInteger("215872848293952797667059038254236574542"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("325629922445073537381874007153931917509"), 128, false),
               BigFixedPoint(BigInteger("98778757057029163154891549982764432298"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, false),
               BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("160407997365957196655387116530659854223"), 128, false),
               BigFixedPoint(BigInteger("300102255270364911803307903244229695504"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("215872848293952797667059038254236574542"), 128, false),
               BigFixedPoint(BigInteger("263041826724899848421390153303749767787"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("98778757057029163154891549982764432298"), 128, true),
               BigFixedPoint(BigInteger("325629922445073537381874007153931917509"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, true),
               BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("300102255270364911803307903244229695504"), 128, true),
               BigFixedPoint(BigInteger("160407997365957196655387116530659854223"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("263041826724899848421390153303749767787"), 128, false),
               BigFixedPoint(BigInteger("215872848293952797667059038254236574542"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("325629922445073537381874007153931917509"), 128, true),
               BigFixedPoint(BigInteger("98778757057029163154891549982764432298"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, true),
               BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("160407997365957196655387116530659854223"), 128, true),
               BigFixedPoint(BigInteger("300102255270364911803307903244229695504"), 128, false)),
    BigComplex(BigFixedPoint(BigInteger("215872848293952797667059038254236574542"), 128, true),
               BigFixedPoint(BigInteger("263041826724899848421390153303749767787"), 128, true)),
    BigComplex(BigFixedPoint(BigInteger("98778757057029163154891549982764432298"), 128, false),
               BigFixedPoint(BigInteger("325629922445073537381874007153931917509"), 128, false)),
};

const auto r_roots_32_scale = 128;

const auto r_roots       = r_roots_32;
const auto r_roots_scale = r_roots_32_scale;
const auto rN            = 32;

BigCMatrix getRU();

BigCMatrix getRUInverse();

// Multiplies input vector by U matrix
std::vector<BigComplex> multU(const BigCMatrix& U, std::vector<BigFixedPoint> input);

std::vector<BigFixedPoint> multUInverse(const BigCMatrix& UInv, std::vector<BigComplex> input);

// Primitive roots of unity for various m
// Here we define exp(2 pi i / m) for m = 8, 16, ..., 262144
// They are roots for X^{m/2} + 1, and offers m/4 C slots
// The minimal m is 8 because the minimal C slot we have is 2 (sparse packing of 1 complex number)
// If we account for Z slot, minimal m is much larger depending on the choice
const auto R_ROOT_M8      = BigComplex(BigFixedPoint(BigInteger("240615969168004511545033772477625056927"), 128, false),
                                       BigFixedPoint(BigInteger("240615969168004511545033772477625056927"), 128, false));
const auto R_ROOT_M16     = BigComplex(BigFixedPoint(BigInteger("314379914072750776175968446001973125505"), 128, false),
                                       BigFixedPoint(BigInteger("130220424146621615521356287377630227996"), 128, false));
const auto R_ROOT_M32     = BigComplex(BigFixedPoint(BigInteger("333743936656827580393945988193699131951"), 128, false),
                                       BigFixedPoint(BigInteger("66385796539016198537004233668372411924"), 128, false));
const auto R_ROOT_M64     = BigComplex(BigFixedPoint(BigInteger("338643814315582355912937410983851649519"), 128, false),
                                       BigFixedPoint(BigInteger("33353504510164655995148692934415085169"), 128, false));
const auto R_ROOT_M128    = BigComplex(BigFixedPoint(BigInteger("339872481907374595982655373237523399998"), 128, false),
                                       BigFixedPoint(BigInteger("16696864359439569162511425913059001623"), 128, false));
const auto R_ROOT_M256    = BigComplex(BigFixedPoint(BigInteger("340179880234010501256410555240815138610"), 128, false),
                                       BigFixedPoint(BigInteger("8350947328924239869306493724758181497"), 128, false));
const auto R_ROOT_M512    = BigComplex(BigFixedPoint(BigInteger("340256744284537774291616451567166656226"), 128, false),
                                       BigFixedPoint(BigInteger("4175788093625692072829030064408394408"), 128, false));
const auto R_ROOT_M1024   = BigComplex(BigFixedPoint(BigInteger("340275961201545362514682404085958700136"), 128, false),
                                       BigFixedPoint(BigInteger("2087933351568136955647909863034040375"), 128, false));
const auto R_ROOT_M2048   = BigComplex(BigFixedPoint(BigInteger("340280765487321862450779202347741264351"), 128, false),
                                       BigFixedPoint(BigInteger("1043971588913162969962874901364098823"), 128, false));
const auto R_ROOT_M4096   = BigComplex(BigFixedPoint(BigInteger("340281966562298792572160534645360975660"), 128, false),
                                       BigFixedPoint(BigInteger("521986408598802148050493966777353034"), 128, false));
const auto R_ROOT_M8192   = BigComplex(BigFixedPoint(BigInteger("340282266831263825696362301831656808750"), 128, false),
                                       BigFixedPoint(BigInteger("260993281067212527299757014848791288"), 128, false));
const auto R_ROOT_M16384  = BigComplex(BigFixedPoint(BigInteger("340282341898518884018790830072534181986"), 128, false),
                                       BigFixedPoint(BigInteger("130496650129583753759140694874978295"), 128, false));
const auto R_ROOT_M32768  = BigComplex(BigFixedPoint(BigInteger("340282360665333511102050687223426738058"), 128, false),
                                       BigFixedPoint(BigInteger("65248326264289096219790833612199600"), 128, false));
const auto R_ROOT_M65536  = BigComplex(BigFixedPoint(BigInteger("340282365357037221779282487371359009126"), 128, false),
                                       BigFixedPoint(BigInteger("32624163282081701561065576602700777"), 128, false));
const auto R_ROOT_M131072 = BigComplex(BigFixedPoint(BigInteger("340282366529963152817741505908074724901"), 128, false),
                                       BigFixedPoint(BigInteger("16312081659782994994230389606789418"), 128, false));
const auto R_ROOT_M262144 = BigComplex(BigFixedPoint(BigInteger("340282366823194635787928202577525532155"), 128, false),
                                       BigFixedPoint(BigInteger("8156040832234265524836811571533929"), 128, false));

const std::map<uint32_t, BigComplex> R_ROOT_MAP = {
    {8, R_ROOT_M8},         {16, R_ROOT_M16},       {32, R_ROOT_M32},         {64, R_ROOT_M64},
    {128, R_ROOT_M128},     {256, R_ROOT_M256},     {512, R_ROOT_M512},       {1024, R_ROOT_M1024},
    {2048, R_ROOT_M2048},   {4096, R_ROOT_M4096},   {8192, R_ROOT_M8192},     {16384, R_ROOT_M16384},
    {32768, R_ROOT_M32768}, {65536, R_ROOT_M65536}, {131072, R_ROOT_M131072}, {262144, R_ROOT_M262144},
};

}  // namespace lbcrypto

#endif