// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.numerics.switching;

import static java.lang.Boolean.parseBoolean;
import static java.lang.Double.parseDouble;
import static java.lang.String.format;
import static java.util.Arrays.copyOfRange;

import java.util.Arrays;

/**
 * Static class responsible for parsing String arrays into univariate switching functions.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class UnivariateFunctionFactory {

  /** Static only class. */
  private UnivariateFunctionFactory() {
  }

  /**
   * Parse an array of Strings terminating in the description of a univariate switching function.
   *
   * @param toks Array of Strings, terminating with the description of a function.
   * @param offset Index of the first relevant token.
   * @return A parsed univariate switching function.
   */
  public static UnivariateSwitchingFunction parseUSF(String[] toks, int offset) {
    return parseUSF(copyOfRange(toks, offset, toks.length));
  }

  /**
   * Parse an array of Strings describing a univariate switching function.
   *
   * @param toks Descriptive Strings.
   * @return A parsed univariate switching function.
   */
  public static UnivariateSwitchingFunction parseUSF(String[] toks) {
    String selectionString = toks[0].toUpperCase()
        .replaceAll("-", "")
        .replaceAll("_", "")
        .replaceAll(" ", "");
    return switch (selectionString) {
      case "BELL", "BELLCURVE", "BELLCURVESWITCH" -> parseBell(toks);
      case "CONSTANT", "NONE", "FLAT" -> new ConstantSwitch();
      case "LINEARDERIVATIVE" -> new LinearDerivativeSwitch();
      case "MULTIPLICATIVE", "PENTICHERMITE" -> parseMultiplicative(toks);
      case "POWER" -> parsePower(toks);
      case "LINEAR" -> parseSpecificPow(1.0, toks);
      case "QUADRATIC" -> parseSpecificPow(2.0, toks);
      case "CUBIC" -> parseSpecificPow(3.0, toks);
      case "TRIGONOMETRIC", "TRIG", "SINSQUARED" -> parseTrig(toks);
      default -> throw new IllegalArgumentException(
          format(" Could not parse %s as a valid univariate switching function!",
              Arrays.toString(toks)));
    };
  }

  private static BellCurveSwitch parseBell(String[] toks) {
    double midpoint = 0.5;
    if (toks.length > 2) {
      midpoint = parseDouble(toks[1]);
    }
    double width = 1.0;
    if (toks.length > 3) {
      width = parseDouble(toks[2]);
    }
    return new BellCurveSwitch(midpoint, width);
  }

  private static MultiplicativeSwitch parseMultiplicative(String[] toks) {
    if (toks.length > 3) {
      return new MultiplicativeSwitch(parseDouble(toks[2]), parseDouble(toks[1]));
    } else {
      return new MultiplicativeSwitch();
    }
  }

  private static PowerSwitch parsePower(String[] toks) {
    double pow = 1.0;
    if (toks.length > 1) {
      pow = parseDouble(toks[1]);
    }
    double alpha = 1.0;
    if (toks.length > 2) {
      alpha = parseDouble(toks[2]);
    }
    return new PowerSwitch(alpha, pow);
  }

  private static PowerSwitch parseSpecificPow(double pow, String[] toks) {
    double alpha = 1.0;
    if (toks.length > 1) {
      alpha = parseDouble(toks[1]);
    }
    return new PowerSwitch(alpha, pow);
  }

  private static SquaredTrigSwitch parseTrig(String[] toks) {
    boolean trig = false;
    if (toks.length > 1) {
      trig = parseBoolean(toks[1]);
    }
    if (toks.length > 2) {
      return new SquaredTrigSwitch(parseDouble(toks[2]), trig);
    } else {
      return new SquaredTrigSwitch(trig);
    }
  }
}
