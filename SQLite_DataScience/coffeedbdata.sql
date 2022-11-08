/*
 * SQLite INSERT statement examples.
 * Created by Alvin Alexander, http://alvinalexander.com
 * Released under the Creative Commons License.
 */

--
-- coffees
--

INSERT INTO coffees VALUES (null, 'Colombian', 7.99);
INSERT INTO coffees VALUES (null, 'French_Roast', 8.99);
INSERT INTO coffees VALUES (null, 'Espresso', 9.99);
INSERT INTO coffees VALUES (null, 'Colombian_Decaf', 8.99);
INSERT INTO coffees VALUES (null, 'French_Roast_Decaf', 9.99);

--
-- salespeople
--

INSERT INTO salespeople VALUES (null, 'Fred', 'Flinstone', 10.0);
INSERT INTO salespeople VALUES (null, 'Barney', 'Rubble', 10.0);

--
-- customers
--

INSERT INTO customers VALUES (null, 'ACME, INC.', '101 Main Street', 'Anchorage', 'AK', '99501');
INSERT INTO customers VALUES (null, 'FOOBAR', '200 Foo Way', 'Louisville', 'KY', '40207');

--
-- orders
--

INSERT INTO orders VALUES (null, 1, 1);
INSERT INTO orders VALUES (null, 2, 2);

--
-- order_items
--

-- insert one order
INSERT INTO order_items VALUES (null, 1, 1, 5);
INSERT INTO order_items VALUES (null, 1, 2, 8);

-- insert a second order
INSERT INTO order_items VALUES (null, 2, 3, 6);
INSERT INTO order_items VALUES (null, 2, 1, 10);