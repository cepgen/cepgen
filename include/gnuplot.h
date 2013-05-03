#ifndef GNUPLOT_H_
#define GNUPLOT_H_

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>

std::string GetRandomString(int);

typedef enum {
  GP_CLASSIC = 0,
  GP_HISTOGRAM = 1,
  GP_HEATMAP = 2
} plot_type;

/**
 * This object allows to invoke gnuplot, the portable command-line driven
 * graphing utility for Linux, OS/2, MS Windows, OSX, VMS, and many other
 * platforms.
 * @brief Plotting utility used in control plots generation
 */
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
    void SetGrid(bool grid_=true);
    /**
     * @brief Toggles the grid for both the axes
     */
    void SetLogy(bool logy_=true);
    /**
     * @brief Sets the title for the graph
     */
    void SetTitle(std::string);
    /**
     * @brief Sets the name for the graph
     */
    void SetName(std::string);
    /**
     * @brief Sets the caption for the \f$x\f$-axis
     */
    void SetXAxisTitle(std::string);
    /**
     * @brief Sets the caption for the \f$y\f$-axis
     */
    void SetYAxisTitle(std::string);
    /**
     * Simplest way to use this object if you know how to plot things with
     * gnuplot : just create the interpreter, and then feed your commands
     * to generate the graph. 
     * @brief Feeds a command line to the Gnuplot interpreter
     * @param &command The Gnuplot-formatted command line to feed
     */
    void operator <<(const std::string &command);
    void SetHistogram(int,double,double,std::string name_="");
    /**
     * Adds an entry to the histogram in case the current object is set as a
     * GP_HISTOGRAM-type plotter
     * @brief Add an entry to the histogram
     * @param value_ The value to insert into the histogram
     * @return An integer less than or equal 0 if an error has been observed, or
     *  - 1 if the value was within the histogram's range
     *  - 2 if the value was set in the underflow bin
     *  - 3 if the value was set in the overflow bin
     */
    int Fill(double);
    /**
     * Adds an entry to the histogram in case the current object is set as a
     * GP_HISTOGRAM-type plotter
     * @brief Add an entry to the histogram
     * @param value_ The value to insert into the histogram
     * @param weight_ The weight of this entry in the histogram
     * @return An integer less than or equal 0 if an error has been observed, or
     *  - 1 if the value was within the histogram's range
     *  - 2 if the value was set in the underflow bin
     *  - 3 if the value was set in the overflow bin
     */
    int Fill(double,double);
    int Fill2(double,double);
    int Fill2(double,double,double);
    int DrawHistogram();
  protected:
    /** @brief The pipe used to feed the Gnuplot interpreter */
    FILE *_pipe;
  private:
    /**
     * The type of graph involved :
     * - GP_CLASSIC for a simple graph (list of \f$(x,y)\f$ pairs)
     * - GP_HISTOGRAM for a frequency histogram (list of \f$(x, f(x))\f$ pairs)
     * - GP_HISTOGRAM for a 2-dimensional histogram (list of \f$(x, y, f(x, y))\f$
     * triplets)
     * @brief Type of graph to draw
     */
    plot_type _type;
    /** @brief Bin bounds for the histogram */
    double *_histBounds;
    /** @brief Values for each bins in the histogram */
    double *_histValues;
    /** @brief Number of bins in the histogram */
    int _histNum;
    /** @brief Lowest value of the histogram range */
    double _histLow;
    /** @brief Highest value of the histogram range */
    double _histHigh;
    /** @brief Underflow bin for the histogram */
    double _histUnderflow;
    /** @brief Overflow bin for the histogram */
    double _histOverflow;
    /** @brief Is the object complete enough to be drawn ? */
    bool _isPlottable;
    /** @brief Graph title */
    std::string _title;
    /** @brief Graph name */
    std::string _name;
    std::string _outputFile;
    std::string _tmpFile;
};

#endif

