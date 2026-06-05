import ROOT
import plot_summary

def main():
    #for mass in ['15','20','30','40','50','55']:
    for mass in ['30']:
        #for lifetime in ["0", "10", "20", "50", "100", "1000"]:
        for lifetime in ["20"]:
            #for year in ["2016", "2017", "2018","2022", "2023","2024"]:
            for year in ["2022"]:
                plot_summary.run(mass, lifetime, year)

if __name__ == "__main__":
    main()
