#ifndef GNUPLOT_H_
#define GNUPLOT_H_

#include <string>
#include <fstream>
#include <iostream>
#include <cstdio>

class Gnuplot {
  public:
    Gnuplot(std::string outFile_="") ;
    ~Gnuplot();
    /**
     * Sets the file on which the graph has to be produced
     */
    void SetOutputFile(std::string);
    /**
     * @brief Toggles the logarithmic scale for the \f$y\f$-axis
     */
    void SetLogy(bool logy_=true);
    /**
     * @brief Sets the title for the graph
     */
    void SetTitle(std::string);
    /**
     * @brief Sets the caption for the \f$x\f$-axis
     */
    void SetXAxisTitle(std::string);
    /**
     * @brief Sets the caption for the \f$y\f$-axis
     */
    void SetYAxisTitle(std::string);
    /**
     * @brief Feeds a command line to the Gnuplot interpreter
     * @param &command The Gnuplot-formatted command line to feed
     */
    void operator <<(const std::string &command);
    void SetHistogram(int,double,double);
    int Fill(double);
    int DrawHistogram();
  protected:
    /** @brief The pipe used to feed the Gnuplot interpreter */
    FILE *_pipe;
  private:
    bool _isHist;
    double *_histBounds;
    int *_histValues;
    int _histNum;
    double _histLow;
    double _histHigh;
    int _histUnderflow;
    int _histOverflow;
};

#endif

