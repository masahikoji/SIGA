import csv, sys, json
src=sys.argv[1]; outdir=sys.argv[2]; tag=sys.argv[3] if len(sys.argv)>3 else "null"
STRESS=(tag=="stress")
suf="-stress" if STRESS else ""
fsuf="_stress" if STRESS else ""
rows=list(csv.DictReader(open(src)))
f=float
designs=["K2_n200","K2_n1000","K5_n400","K5_n2000"]
dlabel={"K2_n200":"2 factors, $n=200$","K2_n1000":"2 factors, $n=1000$","K5_n400":"5 factors, $n=400$","K5_n2000":"5 factors, $n=2000$"}
dshort={"K2_n200":"$F=2$, $n=200$","K2_n1000":"$F=2$, $n=1000$","K5_n400":"$F=5$, $n=400$","K5_n2000":"$F=5$, $n=2000$"}
objlabel={"superiority":"superiority","noninferiority":"non-inferiority","equivalence":"equivalence"}
rolelabel={"type1":"","type1_lower":" (lower limit)","type1_upper":" (upper limit)"}
def pct(x): return "%.2f"%(100*f(x))
def rng(vals): return "%.2f--%.2f"%(100*min(vals),100*max(vals))
# ---- main-text compact table ----
def signed_max_dev(sub,m):
    devs=[100*f(r[m+"_minus_nominal"]) for r in sub]
    k=max(range(len(devs)), key=lambda i: abs(devs[i]))
    return "$%+.2f$"%devs[k]
L=[]
L.append(r"\begin{table}[!htbp]")
L.append(r"\centering")
L.append((r"\caption{Heterogeneous-effect null study under a strong biased coin ($p_{\rm bc}=.95$, outcome standard deviation $.25$, $\max_s|d_s|=1$ for the continuous outcome): largest deviation (percentage points, with sign) of the rejection probability from the nominal level over the four null values and both scores, the mean estimated variance ratio, and the range of the ratio of the variance of the regenerated statistics to $\widehat V_{{\rm R},n}$.}" if STRESS else r"\caption{Heterogeneous-effect null study: largest deviation (percentage points, with sign) of the rejection probability from the nominal level over the four null values and both scores, the mean estimated variance ratio, and the range of the ratio of the variance of the regenerated statistics to $\widehat V_{{\rm R},n}$.}"))
L.append(r"\label{tab:hetero-null"+suf+"}")
L.append(r"\begin{threeparttable}\begin{singlespace}\small\setlength{\tabcolsep}{4pt}")
L.append(r"\begin{tabular}{@{}llrrrrll@{}}")
L.append(r"\toprule")
L.append(r"Design & Outcome & RT & CRT & SIGA-S & SIGA-R & $\widehat V_{\rm R}/\widehat V_{\rm S}$ & $\widehat V^{*}/\widehat V_{\rm R}$\\")
L.append(r"\midrule")
for d in designs:
    for oc in ["continuous","binary"]:
        sub=[r for r in rows if r["design_label"]==d and r["outcome_type"]==oc]
        if not sub: continue
        vr=[f(r["mean_ratio_vr_model_over_vs"]) for r in sub]; ve=[f(r["mean_ratio_vrt_empirical_over_vr_model"]) for r in sub]
        L.append(" & ".join([dlabel[d],oc,signed_max_dev(sub,"rt"),signed_max_dev(sub,"crt_data_db"),signed_max_dev(sub,"siga_s"),signed_max_dev(sub,"siga_r"),
                             "%.3f"%(sum(vr)/len(vr)),"%.3f--%.3f"%(min(ve),max(ve))])+r"\\")
L.append(r"\bottomrule")
L.append(r"\end{tabular}")
L.append(r"\begin{tablenotes}[flushleft]\small")
L.append(r"\item Nominal levels are 5\% (superiority, two-sided; equivalence, two one-sided tests at 5\%) and 2.5\% (non-inferiority); the Monte Carlo standard error of a rejection probability from 100,000 outer trials is 0.07 points at 5\% and 0.05 points at 2.5\%. RT is the reference test and CRT the corrected test with the data-based, noise-adjusted $\widehat{\bm d}(b)$. $\widehat V^{*}$ is the variance of the 4,999 regenerated statistics in a trial, averaged over trials. The negative binary entries for RT and CRT are the unadjusted superiority scenarios, where ties of the lattice-valued score are counted as exceedances. Complete results are in Supplementary Appendix~J.")
L.append(r"\end{tablenotes}\end{singlespace}\end{threeparttable}")
L.append(r"\end{table}")
open(outdir+"/main_table_hetero_null"+fsuf+".tex","w").write("\n".join(L)+"\n")
# ---- supplementary full table (landscape longtable) ----
S=[]
S.append(r"\clearpage")
S.append(r"\begin{landscape}\begin{singlespace}\begingroup\footnotesize\setlength{\tabcolsep}{2pt}\renewcommand{\arraystretch}{1.05}\setlength{\LTleft}{\fill}\setlength{\LTright}{\fill}")
S.append(r"\begin{longtable}{@{}llllllllll@{}}")
S.append((r"\caption{Rejection probabilities (\%) at the null values under heterogeneous treatment effects and a strong biased coin ($p_{\rm bc}=.95$), with Monte Carlo standard errors in parentheses.}" if STRESS else r"\caption{Rejection probabilities (\%) at the null values under heterogeneous treatment effects, with Monte Carlo standard errors in parentheses.}"))
S.append("\\label{tab:supp-hetero-null"+suf+"}\\\\")
S.append(r"\toprule")
S.append(r"Design & Outcome & Null value & Score & RT & CRT (model) & CRT (data) & CRT (data, adj.) & SIGA-S & SIGA-R\\")
S.append(r"\midrule\endfirsthead")
S.append(r"\multicolumn{10}{c}{\tablename\ \thetable{} -- continued}\\\toprule")
S.append(r"Design & Outcome & Null value & Score & RT & CRT (model) & CRT (data) & CRT (data, adj.) & SIGA-S & SIGA-R\\")
S.append(r"\midrule\endhead")
S.append(r"\midrule\multicolumn{10}{r}{Continued on next page}\\\endfoot")
S.append(r"\bottomrule\endlastfoot")
def cell(r,m): return "%s (%s)"%(pct(r[m]),pct(r[m+"_se"]))
for d in designs:
    for oc in ["continuous","binary"]:
        for r in sorted([r for r in rows if r["design_label"]==d and r["outcome_type"]==oc], key=lambda r:(int(r["scenario_id"]), r["analysis"]!="unadjusted")):
            nv="%s%s, $%s$"%(objlabel[r["objective"]],rolelabel[r["role"]],("%.2f"%f(r["true_effect"])).replace("-","-"))
            S.append(" & ".join([dshort[d],oc,nv,r["analysis"],cell(r,"rt"),cell(r,"crt_model"),cell(r,"crt_data"),cell(r,"crt_data_db"),cell(r,"siga_s"),cell(r,"siga_r")])+r"\\")
S.append(r"\end{longtable}")
S.append(r"\noindent{\small Nominal levels: 5\% for superiority (two-sided) and equivalence (two one-sided tests at 5\%), 2.5\% for non-inferiority. All procedures use the same 100,000 outer trials per scenario. The Gaussian approximations are used without the continuity correction of Appendix~D.}")
S.append(r"\endgroup\end{singlespace}\end{landscape}")
open(outdir+"/supp_table_hetero_null"+fsuf+".tex","w").write("\n".join(S)+"\n")
# ---- supplementary diagnostics table ----
D=[]
D.append(r"\begin{table}[!htbp]\centering")
D.append((r"\caption{Variance diagnostics in the heterogeneous-effect null study under a strong biased coin ($p_{\rm bc}=.95$): means over the 100,000 outer trials at the tested null value, averaged over the four null values.}" if STRESS else r"\caption{Variance diagnostics in the heterogeneous-effect null study: means over the 100,000 outer trials at the tested null value, averaged over the four null values.}"))
D.append(r"\label{tab:supp-hetero-diag"+suf+"}")
D.append(r"\begin{threeparttable}\begin{singlespace}\footnotesize\setlength{\tabcolsep}{3pt}")
D.append(r"\begin{tabular}{@{}lllllllll@{}}\toprule")
D.append(r"Design & Outcome & Score & $\widehat\lambda$ (model) & $\widehat\lambda$ (data) & $\widehat\lambda$ (adj.) & $\widehat V_{\rm R}/\widehat V_{\rm S}$ & $\widehat V^{*}/\widehat V_{\rm S}$ & $\widehat V^{*}/\widehat V_{\rm R}$\\\midrule")
for d in designs:
    for oc in ["continuous","binary"]:
        for a in ["unadjusted","adjusted"]:
            sub=[r for r in rows if r["design_label"]==d and r["outcome_type"]==oc and r["analysis"]==a]
            if not sub: continue
            m=lambda k: sum(f(r[k]) for r in sub)/len(sub)
            D.append(" & ".join([dshort[d],oc,a,"%.4f"%m("mean_lambda_model"),"%.4f"%m("mean_lambda_data"),"%.4f"%m("mean_lambda_data_db"),"%.4f"%m("mean_ratio_vr_model_over_vs"),"%.4f"%m("mean_ratio_vrt_empirical_over_vs"),"%.4f"%m("mean_ratio_vrt_empirical_over_vr_model")])+r"\\")
D.append(r"\bottomrule\end{tabular}")
D.append(r"\begin{tablenotes}[flushleft]\small\item $\widehat V^{*}$ is the sample variance of the 4,999 regenerated statistics of a trial. The floor in the SIGA-R variance was never active.\end{tablenotes}")
D.append(r"\end{singlespace}\end{threeparttable}\end{table}")
open(outdir+"/supp_table_hetero_diag"+fsuf+".tex","w").write("\n".join(D)+"\n")
print("tables written")
