// =====================================================================
// IF_Schwellenwert_Hilfe.ijm
// Fiji / ImageJ 1.54t
//
// Hilfsmakro zum Festlegen des Schwellenwerts (ohne Negativkontrolle).
// Fuer jedes Bild werden in der Auswertflaeche (Metastase minus Nekrose /
// Ausschluss, falls ROIs vorhanden) mehrere Schwellenwert-Methoden
// ausprobiert und als Kachelbild (Montage) gespeichert:
//   Kachel 1 = Originalbild (Kontrast gestreckt)
//   weitere Kacheln = Signalmaske je Methode (weiss = Signal), plus optional
//   deine eigenen Kandidaten-Schwellenwerte ("Fest=...").
// Zusaetzlich entsteht eine Tabelle mit dem Schwellenwert jeder Methode pro Bild
// und eine Zusammenfassung (Median/Min/Max pro Methode ueber alle Bilder).
//
// Vorgehen: Ein paar repraesentative Bilder beider Gruppen in einen Ordner legen,
// Makro laufen lassen, Montagen mit dem Originalbild vergleichen, einen festen Wert
// waehlen (z.B. Median der Methode, die Signal und Artefakte am besten trennt),
// ihn als Kandidat nochmal testen und dann denselben Wert fuer ALLE Bilder verwenden.
// Idealerweise verblindet gegenueber der Gruppe.
//
// ROIs: Existiert <Bildname>_ROIs.zip im Ordner Auswertung_IF_Neu, wird es
// verwendet. Sonst zeichnest du einen groben Bereich (wird NICHT gespeichert).
// =====================================================================

requires("1.54t");
setOption("BlackBackground", true);


var useComma = false;
var rgArea = 0;
var rgPosArea = 0;
var rgPct = 0;
var rgRaw = 0;
var rgIntDen = 0;
var rgMeanPos = 0;
var rgMeanAll = 0;
var rgIntDenPerArea = 0;

// ---------- Hilfsfunktionen ----------
function fmt(x, dec) {
    s = d2s(x, dec);
    if (useComma) s = replace(s, "\\.", ",");
    return s;
}

function closeID(id) {
    selectImage(id);
    close();
}

function blankMask(w, h) {
    newImage("blank", "8-bit black", w, h, 1);
    return getImageID();
}

// 8-bit-Maske (255 = innerhalb der ROI) fuer ROI mit Index roiIdx
function makeRoiMask(roiIdx, w, h) {
    newImage("roiMask", "8-bit black", w, h, 1);
    id = getImageID();
    roiManager("select", roiIdx);
    setColor(255);
    fill();
    run("Select None");
    return id;
}

// a OR b (neues Bild); a und b werden geschlossen
function orMask(aID, bID) {
    imageCalculator("Max create", aID, bID);
    newID = getImageID();
    closeID(aID);
    closeID(bID);
    return newID;
}

// a AND b (neues Bild), a und b bleiben offen
function andMask(aID, bID) {
    imageCalculator("AND create", aID, bID);
    return getImageID();
}

// a minus b (neues Bild), a und b bleiben offen
function subMask(aID, bID) {
    imageCalculator("Subtract create", aID, bID);
    return getImageID();
}

function countWhite(id) {
    selectImage(id);
    run("Select None");
    getHistogram(values, counts, 256);
    return counts[255];
}

// Summe der Intensitaeten (32-bit-Kanal chFID) ueber alle weissen Pixel der Maske
function sumMasked(chFID, maskID) {
    selectImage(maskID);
    run("Select None");
    run("Duplicate...", "title=m01");
    m01 = getImageID();
    run("Divide...", "value=255");
    imageCalculator("Multiply create 32-bit", chFID, m01);
    prodID = getImageID();
    getRawStatistics(nPix, meanVal);
    total = meanVal * nPix;
    closeID(prodID);
    closeID(m01);
    return total;
}

// Entfernt Objekte < minPx Pixel aus einer Maske; gibt IMMER ein neues Bild zurueck
function filterSmall(maskID, minPx) {
    selectImage(maskID);
    run("Select None");
    if (minPx <= 1) {
        run("Duplicate...", "title=filt");
        return getImageID();
    }
    run("Analyze Particles...", "size=" + minPx + "-Infinity pixel show=[Count Masks]");
    setThreshold(1, 1e9);
    run("Convert to Mask");
    return getImageID();
}

// Maske der Pixel mit Distanzwert in [loPx, hiPx) (hiPx < 0 = unbegrenzt)
function rangeMask(distID, loPx, hiPx) {
    selectImage(distID);
    run("Select None");
    run("Duplicate...", "title=rng");
    if (hiPx < 0) {
        setThreshold(loPx, 1e9);
    } else {
        setThreshold(loPx, hiPx - 0.0001);
    }
    run("Convert to Mask");
    return getImageID();
}

// Flaechen-/Intensitaetskennzahlen fuer eine Region; Ergebnisse in rg*-Variablen
function measureRegion(regionID, posID, chFID, pw, ph) {
    regPx = countWhite(regionID);
    rgArea = regPx * pw * ph;
    rgPosArea = 0;
    rgPct = 0;
    rgRaw = 0;
    rgIntDen = 0;
    rgMeanPos = 0;
    rgMeanAll = 0;
    rgIntDenPerArea = 0;
    if (regPx > 0) {
        rgMeanAll = sumMasked(chFID, regionID) / regPx;
        posInID = andMask(posID, regionID);
        posPx = countWhite(posInID);
        if (posPx > 0) {
            rgRaw = sumMasked(chFID, posInID);
            rgMeanPos = rgRaw / posPx;
        }
        closeID(posInID);
        rgPosArea = posPx * pw * ph;
        rgPct = 100 * posPx / regPx;
        rgIntDen = rgRaw * pw * ph;
        rgIntDenPerArea = rgIntDen / rgArea;
    }
}

// Randkorrigierte 32-bit-Distanzkarte (Pixel-Einheiten) zum Metastasenrand
function buildDistanceMap(roiIdx, w, h, padPx) {

    newImage("maskTmp", "8-bit black", w, h, 1);
    roiManager("select", roiIdx);
    setColor(255);
    fill();
    run("Select None");

    makeRectangle(0, 0, 1, h);
    run("Duplicate...", "title=leftCol");
    run("Size...", "width=" + padPx + " height=" + h + " depth=1 interpolation=None");

    selectWindow("maskTmp");
    makeRectangle(w - 1, 0, 1, h);
    run("Duplicate...", "title=rightCol");
    run("Size...", "width=" + padPx + " height=" + h + " depth=1 interpolation=None");

    newImage("xPadded", "8-bit black", w + 2 * padPx, h, 1);
    selectWindow("leftCol"); run("Select All"); run("Copy");
    selectWindow("xPadded"); makeRectangle(0, 0, padPx, h); run("Paste");
    selectWindow("maskTmp"); run("Select All"); run("Copy");
    selectWindow("xPadded"); makeRectangle(padPx, 0, w, h); run("Paste");
    selectWindow("rightCol"); run("Select All"); run("Copy");
    selectWindow("xPadded"); makeRectangle(padPx + w, 0, padPx, h); run("Paste");
    run("Select None");
    close("leftCol");
    close("rightCol");
    close("maskTmp");

    wx = w + 2 * padPx;
    selectWindow("xPadded");
    makeRectangle(0, 0, wx, 1);
    run("Duplicate...", "title=topRow");
    run("Size...", "width=" + wx + " height=" + padPx + " depth=1 interpolation=None");

    selectWindow("xPadded");
    makeRectangle(0, h - 1, wx, 1);
    run("Duplicate...", "title=botRow");
    run("Size...", "width=" + wx + " height=" + padPx + " depth=1 interpolation=None");

    newImage("paddedMask", "8-bit black", wx, h + 2 * padPx, 1);
    selectWindow("topRow"); run("Select All"); run("Copy");
    selectWindow("paddedMask"); makeRectangle(0, 0, wx, padPx); run("Paste");
    selectWindow("xPadded"); run("Select All"); run("Copy");
    selectWindow("paddedMask"); makeRectangle(0, padPx, wx, h); run("Paste");
    selectWindow("botRow"); run("Select All"); run("Copy");
    selectWindow("paddedMask"); makeRectangle(0, padPx + h, wx, padPx); run("Paste");
    run("Select None");
    close("topRow");
    close("botRow");
    close("xPadded");

    selectWindow("paddedMask");
    maskPadID = getImageID();
    run("Options...", "iterations=1 count=1 black edm=32-bit");
    run("Distance Map");
    edmID = getImageID();

    makeRectangle(padPx, padPx, w, h);
    run("Duplicate...", "title=distCrop");
    distID = getImageID();

    selectImage(edmID);
    close();
    if (edmID != maskPadID) {
        selectImage(maskPadID);
        close();
    }
    return distID;
}

// 32-bit-Distanzkarte (Pixel) jedes Pixels zum naechsten Pixel der Maske targetID (255 = Zielobjekt)
function distanceToObject(targetID, w, h) {
    selectImage(targetID);
    run("Select None");
    run("Duplicate...", "title=invTarget");
    invID = getImageID();
    run("Invert");
    run("Options...", "iterations=1 count=1 black edm=32-bit");
    run("Distance Map");
    edmID = getImageID();
    distID = edmID;
    if (edmID != invID) {
        selectImage(invID);
        close();
    }
    return distID;
}

// Kontrollbild fuer eine Metastase in das Gesamt-Kontrollbild eintragen
function addToCtrl(ctrlID, metaID, exclID, posID) {
    selectImage(metaID);
    run("Select None");
    run("Duplicate...", "title=invMeta");
    invID = getImageID();
    run("Invert");
    imageCalculator("Min", ctrlID, invID);
    closeID(invID);

    exclInMeta = andMask(exclID, metaID);
    selectImage(exclInMeta);
    run("Multiply...", "value=0.2");
    imageCalculator("Max", ctrlID, exclInMeta);
    closeID(exclInMeta);

    imageCalculator("Max", ctrlID, posID);
}

// ROIs interaktiv einzeichnen (Metastase, Nekrose, Ausschluss)
function drawRois(origID, imgName) {
    metaCount = 0;
    continueDrawing = true;
    while (continueDrawing) {
        selectImage(origID);
        run("Select None");
        setTool("freehand");
        waitForUser("Metastase " + (metaCount + 1) + " einzeichnen: " + imgName,
            "Metastase mit dem Freihand-Werkzeug umranden, dann OK.\n" +
            "(NICHT mit T bestaetigen - das macht das Makro.)\n\n" +
            "Keine (weitere) auswertbare Metastase -> ohne Auswahl auf OK klicken.");
        if (selectionType() == -1) {
            continueDrawing = false;
        } else {
            roiManager("add");
            idx = roiManager("count") - 1;
            roiManager("select", idx);
            roiManager("rename", "Meta_" + (metaCount + 1));
            selectImage(origID);
            run("Select None");

            if (getBoolean("Enthaelt Metastase " + (metaCount + 1) + " einen nekrotischen Kern, den du einzeichnen moechtest?")) {
                selectImage(origID);
                setTool("freehand");
                waitForUser("Nekrose einzeichnen: Metastase " + (metaCount + 1),
                    "Nekrotischen Kernbereich umranden, dann OK.");
                if (selectionType() != -1) {
                    roiManager("add");
                    idxN = roiManager("count") - 1;
                    roiManager("select", idxN);
                    roiManager("rename", "Nekrose_" + (metaCount + 1));
                    selectImage(origID);
                    run("Select None");
                }
            }

            if (getBoolean("Artefakt-Ausschlussbereiche fuer Metastase " + (metaCount + 1) + " einzeichnen (z.B. einzelne stark gefaerbte Zellen)?")) {
                exclCount = 0;
                more = true;
                while (more) {
                    selectImage(origID);
                    setTool("freehand");
                    waitForUser("Ausschlussbereich " + (exclCount + 1) + " einzeichnen",
                        "Bereich umranden, der von der Quantifizierung ausgeschlossen werden soll, dann OK.");
                    if (selectionType() != -1) {
                        roiManager("add");
                        idxE = roiManager("count") - 1;
                        roiManager("select", idxE);
                        roiManager("rename", "Ausschluss_" + (metaCount + 1) + "_" + (exclCount + 1));
                        exclCount++;
                        selectImage(origID);
                        run("Select None");
                    }
                    more = getBoolean("Weiteren Ausschlussbereich einzeichnen?");
                }
            }

            metaCount++;
            continueDrawing = getBoolean("Weitere Metastase in diesem Bild einzeichnen?");
        }
    }
}

// ---------- Einstellungen ----------
inDir  = getDirectory("Ordner mit den .nd2-Dateien waehlen");
outDir = inDir + "Auswertung_IF_Neu" + File.separator;
helpDir = outDir + "Schwellenwert_Hilfe" + File.separator;
File.makeDirectory(outDir);
File.makeDirectory(helpDir);

Dialog.create("Schwellenwert-Hilfe");
Dialog.addString("Marker (fuer Dateinamen):", "CTargetName");
Dialog.addNumber("Kanal Protein of Interest:", 2);
Dialog.addString("Eigene Kandidaten-Schwellenwerte (kommagetrennt, max. 3, leer = keine):", "");
Dialog.addNumber("Mindestgroesse Signalobjekte (um2, 0 = aus):", 0);
Dialog.addNumber("Gauss-Weichzeichnung nur fuer Maske (sigma in Pixel, 0 = aus):", 0);
Dialog.addNumber("Nur jedes n-te Bild verwenden (1 = alle):", 1);
Dialog.addCheckbox("16-bit-Histogramm nutzen (bei Fehlermeldung abwaehlen)", true);
Dialog.addCheckbox("Dezimaltrennzeichen Komma statt Punkt", false);
Dialog.show();

marker   = Dialog.getString();
candText = Dialog.getString();
pocCh    = Dialog.getNumber();
minUm2   = Dialog.getNumber();
sigma    = Dialog.getNumber();
stepN    = Dialog.getNumber();
use16    = Dialog.getCheckbox();
useComma = Dialog.getCheckbox();
if (stepN < 1) stepN = 1;

autoOpt = " dark";
if (use16) autoOpt = " dark 16-bit";

methods = newArray("Otsu", "Triangle", "Li", "MaxEntropy", "Yen", "Moments", "Huang", "RenyiEntropy");

cands = newArray(0);
if (lengthOf(candText) > 0) {
    parts = split(candText, ",");
    for (i = 0; i < parts.length; i++) {
        v = parseFloat(parts[i]);
        if (!isNaN(v) && cands.length < 3) cands = Array.concat(cands, v);
    }
}

allFiles = getFileList(inDir);
files = newArray(0);
for (i = 0; i < allFiles.length; i++) {
    if (endsWith(toLowerCase(allFiles[i]), ".nd2"))
        files = Array.concat(files, allFiles[i]);
}
if (files.length == 0) exit("Keine .nd2-Dateien im gewaehlten Ordner gefunden.");

nImgMax = files.length;
thrAll = newArray(methods.length * nImgMax);
procCount = 0;

outTab = "Bild\tMethode\tSchwellenwert\tProzent_Flaeche_positiv_in_Auswertflaeche\tMindestgroesse_um2\n";

run("Close All");
roiManager("reset");
run("ROI Manager...");

for (f = 0; f < files.length; f = f + stepN) {
    name = files[f];
    baseName = substring(name, 0, lastIndexOf(name, "."));
    print("(" + (f + 1) + "/" + files.length + ") " + name);

    run("Bio-Formats Importer", "open=[" + inDir + name + "] color_mode=Composite view=Hyperstack stack_order=XYCZT");
    origID = getImageID();
    getDimensions(w, h, nC, nZ, nT);
    getPixelSize(pxUnit, pixelW, pixelH, pixelD);

    if (nC < pocCh) {
        print("   -> uebersprungen: nur " + nC + " Kanaele.");
        run("Close All");
        continue;
    }
    if (nZ > 1) {
        run("Z Project...", "projection=[Max Intensity]");
        zID = getImageID();
        selectImage(origID);
        close();
        origID = zID;
    }
    if (nC > 1) Stack.setDisplayMode("composite");

    // ---------- Auswertbereich ----------
    roiPath = outDir + baseName + "_ROIs.zip";
    roiManager("reset");
    if (File.exists(roiPath)) {
        roiManager("open", roiPath);
        print("   -> vorhandene ROIs verwendet.");
    } else {
        selectImage(origID);
        run("Select None");
        setTool("freehand");
        waitForUser("Bereich einzeichnen: " + name,
            "Metastase (oder groben repraesentativen Bereich) umranden, dann OK.\n" +
            "Dieser Bereich wird NICHT gespeichert.");
        if (selectionType() != -1) {
            roiManager("add");
            roiManager("select", roiManager("count") - 1);
            roiManager("rename", "Meta_1");
            selectImage(origID);
            run("Select None");
        }
    }
    nTotal = roiManager("count");
    if (nTotal == 0) {
        print("   -> kein Bereich, uebersprungen.");
        run("Close All");
        continue;
    }

    // Vereinigung aller Metastasen minus Nekrose/Ausschluss
    metaU = blankMask(w, h);
    exclU = blankMask(w, h);
    for (i = 0; i < nTotal; i++) {
        selectImage(origID);
        roiManager("select", i);
        rname = Roi.getName();
        run("Select None");
        if (startsWith(rname, "Meta_")) {
            mk = makeRoiMask(i, w, h);
            metaU = orMask(metaU, mk);
        } else if (startsWith(rname, "Nekrose_") || startsWith(rname, "Ausschluss_")) {
            ek = makeRoiMask(i, w, h);
            exclU = orMask(exclU, ek);
        }
    }
    validID = subMask(metaU, exclU);
    validPx = countWhite(validID);
    if (validPx == 0) {
        print("   -> leere Auswertflaeche, uebersprungen.");
        run("Close All");
        continue;
    }

    selectImage(origID);
    run("Select None");
    run("Duplicate...", "title=POI duplicate channels=" + pocCh);
    chID = getImageID();

    minPx = round(minUm2 / (pixelW * pixelH));
    sc = 700 / maxOf(w, h);

    nTiles = 1 + methods.length + cands.length;
    tileIDs = newArray(nTiles);
    tileLabels = newArray(nTiles);

    // --- Kachel 0: Original ---
    selectImage(chID);
    run("Select None");
    run("Duplicate...", "title=origTile");
    run("Enhance Contrast", "saturated=0.35");
    run("8-bit");
    srcTile = getImageID();
    run("Scale...", "x=" + sc + " y=" + sc + " interpolation=Bilinear average create");
    tileIDs[0] = getImageID();
    closeID(srcTile);
    tileLabels[0] = "Original";

    // --- Methoden und Kandidaten ---
    for (t = 0; t < methods.length + cands.length; t++) {
        if (t < methods.length) {
            selectImage(validID);
            run("Create Selection");
            selectImage(chID);
            run("Restore Selection");
            setAutoThreshold(methods[t] + autoOpt);
            getThreshold(thrVal, dummyUpper);
            resetThreshold();
            run("Select None");
            lbl = methods[t];
        } else {
            thrVal = cands[t - methods.length];
            lbl = "Fest";
        }

        selectImage(chID);
        run("Select None");
        run("Duplicate...", "title=src");
        if (sigma > 0) run("Gaussian Blur...", "sigma=" + sigma);
        setThreshold(thrVal, 1e9);
        run("Convert to Mask");
        sigRaw = getImageID();
        sigValid = andMask(sigRaw, validID);
        closeID(sigRaw);
        posID = filterSmall(sigValid, minPx);
        closeID(sigValid);

        pct = 100 * countWhite(posID) / validPx;

        newImage("ctrlTile", "8-bit black", w, h, 1);
        ctrlTile = getImageID();
        run("Add...", "value=90");
        addToCtrl(ctrlTile, metaU, exclU, posID);
        closeID(posID);

        selectImage(ctrlTile);
        run("Scale...", "x=" + sc + " y=" + sc + " interpolation=None create");
        tileIDs[1 + t] = getImageID();
        closeID(ctrlTile);

        tileLabels[1 + t] = lbl + " T=" + d2s(thrVal, 0) + " (" + d2s(pct, 1) + "%)";

        if (t < methods.length) {
            thrAll[t * nImgMax + procCount] = thrVal;
            outTab = outTab + name + "\t" + methods[t] + "\t" + fmt(thrVal, 1) + "\t" + fmt(pct, 3) + "\t" + fmt(minUm2, 3) + "\n";
        } else {
            outTab = outTab + name + "\tFest\t" + fmt(thrVal, 1) + "\t" + fmt(pct, 3) + "\t" + fmt(minUm2, 3) + "\n";
        }
    }

    // --- Montage ---
    selectImage(tileIDs[0]);
    tileW = getWidth();
    tileH = getHeight();
    cols = 4;
    rows = floor((nTiles + cols - 1) / cols);
    newImage("montage", "8-bit black", cols * tileW, rows * tileH, 1);
    montID = getImageID();
    fs = 22;
    for (i = 0; i < nTiles; i++) {
        px = (i % cols) * tileW;
        py = floor(i / cols) * tileH;
        selectImage(tileIDs[i]);
        run("Select All");
        run("Copy");
        selectImage(montID);
        makeRectangle(px, py, tileW, tileH);
        run("Paste");
        run("Select None");
        setFont("SansSerif", fs, "bold");
        tw = getStringWidth(tileLabels[i]);
        setColor(0);
        fillRect(px + 4, py + 4, tw + 16, fs + 14);
        setColor(255);
        drawString(tileLabels[i], px + 12, py + fs + 8);
    }
    selectImage(montID);
    saveAs("PNG", helpDir + "Montage_" + baseName + "_" + marker + ".png");

    procCount++;
    File.saveString(outTab, helpDir + "Schwellenwerte_" + marker + ".txt");

    run("Clear Results");
    run("Close All");
    roiManager("reset");
}

// ---------- Zusammenfassung ueber alle Bilder ----------
outSum = "Methode\tAnzahl_Bilder\tMedian\tMinimum\tMaximum\n";
for (t = 0; t < methods.length; t++) {
    if (procCount > 0) {
        vals = newArray(procCount);
        for (i = 0; i < procCount; i++) vals[i] = thrAll[t * nImgMax + i];
        Array.getStatistics(vals, vmin, vmax, vmean, vstd);
        sorted = Array.copy(vals);
        Array.sort(sorted);
        if (procCount % 2 == 1) {
            med = sorted[(procCount - 1) / 2];
        } else {
            med = (sorted[procCount / 2 - 1] + sorted[procCount / 2]) / 2;
        }
        outSum = outSum + methods[t] + "\t" + procCount + "\t" + fmt(med, 1) + "\t" + fmt(vmin, 1) + "\t" + fmt(vmax, 1) + "\n";
        print(methods[t] + ": Median " + d2s(med, 1) + " (Min " + d2s(vmin, 1) + ", Max " + d2s(vmax, 1) + ")");
    }
}
File.saveString(outSum, helpDir + "Schwellenwert_Zusammenfassung_" + marker + ".txt");

print("Fertig. Montagen und Tabellen in: " + helpDir);
showMessage("Fertig", "Schwellenwert-Hilfe abgeschlossen.\n\nMontagen und Tabellen liegen in:\n" + helpDir);
